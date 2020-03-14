 /*********************************************************
 *
 * This source code is part of the Carnegie Mellon Robot
 * Navigation Toolkit (CARMEN)
 *
 * CARMEN Copyright (c) 2002 Michael Montemerlo, Nicholas
 * Roy, Sebastian Thrun, Dirk Haehnel, Cyrill Stachniss,
 * and Jared Glover
 *
 * CARMEN is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation;
 * either version 2 of the License, or (at your option)
 * any later version.
 *
 * CARMEN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General
 * Public License along with CARMEN; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place,
 * Suite 330, Boston, MA  02111-1307 USA
 *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
/* #include <sys/time.h> */
/* #include <unistd.h> */
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "vasco.h"
#include "vascocore_intern.h"

/* 指定された格子の占有確率を返す */
double
get_mapval( int pos_x, int pos_y, carmen_vascocore_map_t map )
{
  if ( pos_x>=0 && pos_x<map.mapsize.x &&
       pos_y>=0 && pos_y<map.mapsize.y ) {
    /* 微小な正の値を足してから返す */
    return( EPSILON + map.mapprob[pos_x][pos_y] );
  } else {
    /* インデックスが範囲外であれば, デフォルトの値を返す */
    return( EPSILON + carmen_vascocore_settings.local_map_std_val );
  }
  return(0);
}

/* 指定された格子のヒットとミスの合計回数を返す */
int
get_mapsum( int pos_x, int pos_y, carmen_vascocore_map_t map )
{
  if ( pos_x>=0 && pos_x<map.mapsize.x &&
       pos_y>=0 && pos_y<map.mapsize.y ) {
    return( map.mapsum[pos_x][pos_y] );
  } else {
    return( 0 );
  }
  return(0);
}

/* 指定された格子のヒット回数を返す */
float
get_maphit( int pos_x, int pos_y, carmen_vascocore_map_t map )
{
  if ( pos_x>=0 && pos_x<map.mapsize.x &&
       pos_y>=0 && pos_y<map.mapsize.y ) {
    return( map.maphit[pos_x][pos_y] );
  } else {
    return( 0 );
  }
  return(0);
}

int
max( int a, int b )
{
  return( b>a?b:a );
}

/* centerとインデックス(x, y)の相対的な位置関係を返す */
int
find_quadrant( carmen_svec2_t center, int x, int y )
{
  /* 右上: 0, 左上: 1, 左下: 2, 右下: 3 */
  if (x>(center.x)) {
    if (y>(center.y)) {
      return(0);
    } else {
      return(3);
    }
  } else {
    if (y>(center.y)) {
      return(1);
    } else {
      return(2);
    }
  }
}

/* centerとの相対的な位置関係を基に,
 * centerからstepsize分だけずらしたインデックスを返す */
carmen_svec2_t
newcenter( carmen_svec2_t center, int i, short stepsize )
{
  /* 右上: 0, 左上: 1, 左下: 2, 右下: 3 */
  carmen_svec2_t ncenter = center;
  switch(i) {
  case 0:
    ncenter.x += stepsize;
    ncenter.y += stepsize;
    break;
  case 1:
    ncenter.x -= stepsize;
    ncenter.y += stepsize;
    break;
  case 2:
    ncenter.x -= stepsize;
    ncenter.y -= stepsize;
    break;
  case 3:
    ncenter.x += stepsize;
    ncenter.y -= stepsize;
    break;
  }
  return(ncenter);
}

/* 4分木のメモリを確保
 * treeは4分木のルートノード, levelは4分木の深さ,
 * centerはルートノードの中心のインデックス,
 * stepsizeは4分木がカバーする地図の領域の半径(格子数)
 * 2のlevel乗がstepsizeに相当 */
void
alloc_tree( carmen_vascocore_quad_tree_t * tree,
	    int level, carmen_svec2_t center, short stepsize )
{
  int i;
  short nstepsize = stepsize/2;
  /* 現在のノードの中心のインデックス, 深さ, 使用状態を設定 */
  tree->center = center;
  tree->level  = level;
  tree->inuse  = FALSE;
  /* 親ノードから子ノードへと再帰的に初期化 */
  if (level>0) {
    /* 右上, 左上, 左下, 右下の子ノードを順に初期化 */
    for( i=0; i<4; i++) {
      tree->elem[i] = (carmen_vascocore_quad_tree_t *)
	      malloc( sizeof(carmen_vascocore_quad_tree_t) );
      carmen_test_alloc(tree->elem[i]);
      /* 各ノードの子ノードを再帰的に初期化
       * 子ノードがカバーする領域の一辺のサイズは, 親ノードの半分 */
      /* newcenter()関数により, 各位置に対応するインデックスを計算 */
      /* level - 1が0, nstepsizeが1となった段階で終了 */
      alloc_tree( tree->elem[i], level-1,
		    newcenter( center, i, nstepsize ), nstepsize );
    }
  }
}

/* 4分木を初期化
 * size_x, size_yは, 4分木がカバーする地図の領域のサイズ(格子数) */
void
initialize_qtree( carmen_vascocore_quad_tree_t * tree, int size_x, int size_y)
{
  /* 4分木の深さは, log2(size_x)またはlog2(size_y)とする */
  int i,v,nlevel = max( (int) ceil(log10(size_x)/log10(2)),
			(int) ceil(log10(size_y)/log10(2)) );
  carmen_svec2_t center;
  if (carmen_vascocore_settings.verbose)
    fprintf( stderr, "* INFO: num levels       = %d\n", nlevel );
  /* vは2の累乗であり, max(size_x, size_y)よりも大きな最小値
   * 例えば, size_xが12, size_yが23であれば, vは32となる
   * このとき, 32個四方の格子をカバーする4分木が作成される */
  v = 1;
  for (i=0;i<nlevel;i++) v=v*2;
  if (carmen_vascocore_settings.verbose) {
    fprintf( stderr, "* INFO: size             = %d/%d\n", size_x, size_y );
    fprintf( stderr, "* INFO: poss. max size   = %d/%d\n", v, v );
  }
  /* ルートノードの中心のインデックスを2倍した値を設定 */
  center.x = v-1;
  center.y = v-1;
  if (carmen_vascocore_settings.verbose) {
    fprintf( stderr, "* INFO: tree center:       %5.1f %5.1f\n",
	     center.x/2.0, center.y/2.0 );
    fprintf( stderr, "* INFO: tree step:         %5.1f %5.1f\n",
	     0.5*(v/2), 0.5*(v/2) );
    fprintf( stderr, "* INFO: allocate tree: ... " );
  }
  /* ルートノードを基点として, 4分木のメモリを確保 */
  alloc_tree( tree, nlevel, center, v );
  if (carmen_vascocore_settings.verbose)
    fprintf( stderr, "done\n" );
}

/* 地図を初期化
 * sx, syは地図の半径(横方向と縦方向の格子数の半分)
 * center_x, center_yは地図の中心の座標
 * resolutionは各格子のサイズ */
void
initialize_map( carmen_vascocore_map_t * map,
		int sx, int sy, int center_x, int center_y,
		double resolution, carmen_point_t start )
{
  int x, y;

  /* 地図のサイズ(格子数)を設定 */
  map->mapsize.x  = sx;
  map->mapsize.y  = sy;
  /* 地図の各格子のサイズを設定 */
  map->resolution = resolution;
  /* TODO */
  map->offset     = start;
  if (carmen_vascocore_settings.verbose) {
    fprintf( stderr, "* INFO: allocating memory ... " );
  }

  /* 地図の各メンバのメモリを確保 */
  map->updated  = carmen_mdalloc( 2, sizeof(unsigned char),  sx, sy );
  carmen_test_alloc(map->updated);
  map->maphit   = carmen_mdalloc( 2, sizeof(float),  sx, sy );
  carmen_test_alloc(map->maphit);
  map->mapsum   = carmen_mdalloc( 2, sizeof(short),  sx, sy );
  carmen_test_alloc(map->mapsum);
  map->mapprob  = carmen_mdalloc( 2, sizeof(float), sx, sy );
  carmen_test_alloc(map->mapprob);
  map->calc     = carmen_mdalloc( 2, sizeof(float), sx, sy );
  carmen_test_alloc(map->calc);
  if (carmen_vascocore_settings.verbose) {
    fprintf( stderr, "done\n" );
  }

  /* 地図の中心の座標を設定 */
  map->center.x = center_x;
  map->center.y = center_y;

  if (carmen_vascocore_settings.verbose) {
    fprintf( stderr, "* INFO: map:            %d %d\n",
	     map->mapsize.x, map->mapsize.y );
    fprintf( stderr, "* INFO: center:         %.1f %.1f\n",
	     map->center.x, map->center.y );
    fprintf( stderr, "* INFO: resolution:      %.2f\n",
	     map->resolution );
    fprintf( stderr, "* INFO: real-size:      [%.1f %.1f] [%.1f %.1f]\n",
	     -sx*map->resolution, sx*map->resolution,
	     -sy*map->resolution, sy*map->resolution );
    fprintf( stderr, "***************************************\n" );
  }

  for (x=0;x<sx;x++) {
    for (y=0;y<sy;y++) {
      /* 各格子の占有確率を初期化 */
      map->mapprob[x][y] = carmen_vascocore_settings.local_map_std_val;
      /* 各格子の畳み込みの計算結果を初期化 */
      map->calc[x][y]    = carmen_vascocore_settings.local_map_std_val;
      /* 各格子のヒット回数を初期化 */
      map->maphit[x][y]  = 0.0;
      /* 各格子のヒットとミスの合計回数を初期化 */
      map->mapsum[x][y]  = 0;
      /* 各格子の畳み込みの計算状態(calcの状態)を初期化
       * 横方向の畳み込みが終了した段階ではUPDT_X,
       * 横方向と縦方向の双方の計算が終了した段階でUPDT_Yとなる */
      map->updated[x][y]  = UPDT_NOT;
    }
  }

  /* 地図全体をカバーする4分木を初期化 */
  initialize_qtree( &(map->qtree), sx, sy );
}

/* 指定された格子の占有確率を計算 */
void
compute_prob_point( carmen_vascocore_map_t *map, int x, int y )
{

  if ( x>=0 && x<map->mapsize.x &&
       y>=0 && y<map->mapsize.y ) {
    /* 各格子に, スキャンデータのヒットまたはミスが1回以上発生した場合 */
    if (get_mapsum(x,y,*map)>0) {
      if (0) {
        /* その格子の占有確率を, ヒットした回数および, 
         * ヒットおよびミスした回数との比率に設定 */
	      map->mapprob[x][y]     = map->maphit[x][y] / map->mapsum[x][y];
      } else {
        /* その格子の占有確率を1に設定 */
	      map->mapprob[x][y]     = 1.0;
      }
    } else {
      /* 各格子に1度も, スキャンデータのヒットまたはミスが起きていない場合
       * その格子の占有確率をデフォルト値に設定 */
      map->mapprob[x][y]     = carmen_vascocore_settings.local_map_std_val;
    }
  }
}

/* 指定された格子を中心とする, ガウスカーネルを適用
 * (x - hk, y)から(x + hk, y)までの横方向の格子に適用
 * std_valは, 未観測である格子のデフォルトの占有確率 */
void
convolve_calc_point( carmen_vascocore_map_t *map,
		     carmen_gauss_kernel_t kernel, int hk,
		     int x, int y, double std_val )
{
  int                     k;
  double                  ksum = 0.0;

  /* xを中心とする半径hkの領域が, 地図の内部に収まる場合 */
  if (x-hk>=0 && x+hk<map->mapsize.x) {
    for (k=0;k<2*hk+1;k++) {
      /* 格子が観測された場合(ヒットまたはミスが発生した場合) */
      if (get_mapsum(x+k-hk,y,*map)>0) {
	      ksum +=  ( kernel.val[k] * map->mapprob[x+k-hk][y] );
      } else {
        /* 格子が観測されていない場合
         * 占有確率としてstd_valで指定されたデフォルト値を使用 */
	      ksum +=  ( kernel.val[k] * std_val );
	    }
    }
    /* 指定された格子に計算結果を格納
     * (x, y)の格子に対するcalcの値は, (x - hk, y)から(x + hk, y)までの
     * 格子にガウスカーネルを適用することで計算される */
    map->calc[x][y]     = ksum;
    /* (x, y)の格子に対するcalcの値は, 横方向の格子を使って計算された */
    map->updated[x][y]  = UPDT_X;
  }
}

/* 指定された格子を中心とする, ガウスカーネルを適用
 * (x, y - hk)から(x, y + hk)までの縦方向の格子に適用
 * hkはガウスカーネルの半径のサイズ(格子数) */
void
convolve_prob_point( carmen_vascocore_map_t *map,
		     carmen_gauss_kernel_t kernel, int hk, int x, int y )
{
  int                     k;
  double                  ksum = 0.0;

  /* yを中心とする半径hkの領域にガウスカーネルを適用 */
  /* y - hk >= 0 && y + hk < map->mapsize.yであることは, 呼び出し側で保証 */
  for (k=0;k<2*hk+1;k++) {
    /* x, yが地図の範囲内であれば, ガウスカーネルを計算 */
    if ( x>=0 && x<map->mapsize.x &&
	       y>=0 && y<map->mapsize.y)
      ksum +=  ( kernel.val[k] * map->calc[x][y+k-hk] );
  }
  /* 指定された格子に計算結果を格納
   * (x, y)の格子に対するmapprobの値は, 横方向にガウスカーネルを適用した後,
   * 縦方向に再び適用することで計算される */
  map->mapprob[x][y]  = ksum;
  /* (x, y)の格子に対するcalcの値は, 横方向かつ縦方向の格子を使って計算された */
  map->updated[x][y]  = UPDT_Y;
}

/* 4分木の使用状態の各ノードに対して,
 * ノードの中心のインデックスに対応する格子の占有確率を計算
 * compute_prob_point()関数を再帰的に呼び出す */
void
compute_prob_treemap( carmen_vascocore_quad_tree_t *tree,
		      carmen_vascocore_map_t *map )
{
  /* 子ノードがいる場合 */
  if ((tree->level)>0 ) {
    /* 使用状態の子ノードに対して, 再帰的に関数を呼び出す */
    if (tree->elem[0]->inuse)
      compute_prob_treemap( tree->elem[0], map );
    if (tree->elem[1]->inuse)
      compute_prob_treemap( tree->elem[1], map );
    if (tree->elem[2]->inuse)
      compute_prob_treemap( tree->elem[2], map );
    if (tree->elem[3]->inuse)
      compute_prob_treemap( tree->elem[3], map );
  } else {
    /* 葉ノードであれば, そのノードがカバーする地図の領域の
     * 中心に位置する格子の占有確率を計算 */
    /* ノードがもつ中心のインデックスは, 実際の値の2倍であることに注意 */
    compute_prob_point( map, (tree->center.x/2), (tree->center.y/2) );
  }
}

/* 4分木の使用状態の各ノードに対して, 横方向にガウスカーネルを適用
 * convolve_calc_point()関数を再帰的に呼び出す */
void
convolve_calc_treemap( carmen_vascocore_quad_tree_t *tree,
		       carmen_vascocore_map_t *map,
		       carmen_gauss_kernel_t kernel, int hk, double std )
{
  int i, j;
  /* 子ノードがいる場合 */
  if ((tree->level)>0 ) {
    /* 使用状態の子ノードに対して, 再帰的に関数を呼び出す */
    if (tree->elem[0]->inuse)
      convolve_calc_treemap( tree->elem[0], map, kernel, hk, std );
    if (tree->elem[1]->inuse)
      convolve_calc_treemap( tree->elem[1], map, kernel, hk, std );
    if (tree->elem[2]->inuse)
      convolve_calc_treemap( tree->elem[2], map, kernel, hk, std );
    if (tree->elem[3]->inuse)
      convolve_calc_treemap( tree->elem[3], map, kernel, hk, std );
  } else {
    /* 葉ノードであれば, その周囲の半径hkの領域に存在する各格子について
     * convolve_calc_point()関数を呼び出して, 横方向にガウスカーネルを適用 */
    /* 葉ノードがカバーする領域の中心のインデックスが,
     * hk + 2以上かつmapsize - hk - 2以下である場合 */
    if ( (tree->center.x/2)>hk+1 && (tree->center.x/2)<map->mapsize.x-hk-1 &&
	       (tree->center.y/2)>hk+1 && (tree->center.y/2)<map->mapsize.y-hk-1 ) {
      for (i=(tree->center.x/2)-hk;i<(tree->center.x/2)+hk;i++) {
        for (j=(tree->center.y/2)-hk;j<(tree->center.y/2)+hk;j++) {
          /* 各格子がまだ, 隣接する横方向の格子によって畳み込まれていない場合 */
          /* 横方向の格子を使ってガウスカーネルを適用 */
          if (map->updated[i][j] != UPDT_X)
            convolve_calc_point( map, kernel, hk, i, j, std );
        }
      }
    }
  }
}

/* 4分木の使用状態の各ノードに対して, 縦方向にガウスカーネルを適用
 * convolve_prob_point()関数を再帰的に呼び出す */
void
convolve_prob_treemap( carmen_vascocore_quad_tree_t *tree,
		       carmen_vascocore_map_t *map,
		       carmen_gauss_kernel_t kernel, int hk )
{
  int i, j;
  /* 子ノードがいる場合 */
  if ((tree->level)>0 ) {
    /* 使用状態の子ノードに対して, 再帰的に関数を呼び出す */
    if (tree->elem[0]->inuse)
      convolve_prob_treemap( tree->elem[0], map, kernel, hk );
    if (tree->elem[1]->inuse)
      convolve_prob_treemap( tree->elem[1], map, kernel, hk );
    if (tree->elem[2]->inuse)
      convolve_prob_treemap( tree->elem[2], map, kernel, hk );
    if (tree->elem[3]->inuse)
      convolve_prob_treemap( tree->elem[3], map, kernel, hk );
  } else {
    /* 葉ノードであれば, その周囲の半径hkの領域に存在する各格子について
     * convolve_prob_point()関数を呼び出して, 縦方向にガウスカーネルを適用 */
    /* 葉ノードがカバーする領域の中心のインデックスが,
     * hk + 2以上かつmapsize - hk - 2以下である場合 */
    if ( (tree->center.x/2)>hk+1 && (tree->center.x/2)<map->mapsize.x-hk-1 &&
	       (tree->center.y/2)>hk+1 && (tree->center.y/2)<map->mapsize.y-hk-1 ) {
      for (i=(tree->center.x/2)-hk;i<(tree->center.x/2)+hk;i++) {
        for (j=(tree->center.y/2)-hk;j<(tree->center.y/2)+hk;j++) {
          /* 各格子がまだ, 隣接する縦方向の格子によって畳み込まれていない場合 */
          /* 縦方向の格子を使ってガウスカーネルを適用 */
          if (map->updated[i][j] != UPDT_Y)
            convolve_prob_point( map, kernel, hk, i, j );
        }
      }
    }
  }
}

/* 4分木の使用状態の各ノードに対して, ガウスカーネルを適用
 * 最初に横方向, 次に縦方向に畳み込みを実行 */
void
convolve_treemap( carmen_vascocore_map_t *map )
{
  int                              i;
  static int                       hk, first_time = TRUE;
  static carmen_gauss_kernel_t     kernel;

  /* 最初の関数呼び出しであれば, ガウスカーネルの各係数を計算 */
  if (first_time) {
    /* ガウスカーネルの半径(local_map_kernel_lenは奇数と想定) */
    hk = (carmen_vascocore_settings.local_map_kernel_len-1)/2;
    /* 1方向のガウスカーネルの係数を計算 */
    kernel = carmen_gauss_kernel( carmen_vascocore_settings.local_map_kernel_len );
    first_time = FALSE;
  }

  /* 各格子の占有確率を計算 */
  compute_prob_treemap(  &(map->qtree), map );
  /* 指定された回数だけ, 畳み込み(ガウスカーネルの適用)を実行 */
  for (i=0;i<carmen_vascocore_settings.local_map_num_convolve;i++) {
    /* 横方向にガウスカーネルを適用 */
    convolve_calc_treemap( &(map->qtree), map, kernel, hk,
			   carmen_vascocore_settings.local_map_std_val );
    /* 縦方向にガウスカーネルを適用 */
    convolve_prob_treemap( &(map->qtree), map, kernel, hk );
  }
}

/* 指定された姿勢と対応する, 地図の格子のインデックスを計算
 * 姿勢の座標系は, 地図の中心を原点とする */
int
compute_map_pos_from_pos( carmen_point_t pos,
			  carmen_vascocore_map_t map, carmen_ivec2_t *v )
{
  v->x = map.center.x + (int) (pos.x/(double)map.resolution);
  v->y = map.center.y + (int) (pos.y/(double)map.resolution);

  /* 地図の範囲内になければ0を返す */
  if (v->x<0 || v->x>map.mapsize.x-1) {
    return(FALSE);
  }
  if (v->y<0 || v->y>map.mapsize.y-1) {
    return(FALSE);
  }

  /* 地図の範囲内であれば1を返す */
  return(TRUE);
}

/* 指定された座標と対応する, 地図の格子のインデックスを計算
 * 座標系は, 地図の中心を原点とする */
int
compute_map_pos_from_vec2( carmen_vec2_t vec,
			   carmen_vascocore_map_t map, carmen_ivec2_t *v )
{
  v->x = map.center.x + (int) (vec.x/(double)map.resolution);
  v->y = map.center.y + (int) (vec.y/(double)map.resolution);

  /* 地図の範囲内になければ0を返す */
  if (v->x<0) {
    return(FALSE);
  } else if (v->x>map.mapsize.x-1) {
    return(FALSE);
  }
  if (v->y<0) {
    return(FALSE);
  } else if (v->y>map.mapsize.y-1) {
    return(FALSE);
  }

  /* 地図の範囲内であれば1を返す */
  return(TRUE);
}

/* 指定した格子に, スキャンデータのヒット回数を加算
 * ヒット回数は浮動小数点数で表現される */
void
mark_maphitpoint( carmen_vascocore_quad_tree_t *tree,
		  carmen_vascocore_map_t *map, int x, int y, float value )
{
  /* 指定された格子を含むノードを, 使用状態に設定 */
  tree->inuse=TRUE;
  /* 子ノードがいる場合は, 再帰的に関数を呼び出す */
  if ((tree->level)>0) {
    /* find_quadrant()関数により, 格子(x, y)を含む適切な子ノードを選択 */
    mark_maphitpoint( tree->elem[find_quadrant( tree->center, x, y )],
		      map, x, y, value  );
  } else {
    /* 葉ノードであればその中心のインデックス,
     * 即ち(x, y)に対応する格子にヒット回数を加算 */
    map->maphit[tree->center.x/2][tree->center.y/2] += value;
  }
}

/* 指定された格子に, スキャンデータの参照回数を記録
 * 格子のmapsumメンバ(ヒットとミスの合計回数)が1だけ加算される */
void
mark_mapsumpoint( carmen_vascocore_quad_tree_t *tree,
		  carmen_vascocore_map_t *map, int x, int y )
{
  /* 指定された格子を含むノードを, 使用状態に設定 */
  tree->inuse=TRUE;
  /* 子ノードがいる場合は, 再帰的に関数を呼び出す */
  if ((tree->level)>0) {
    /* find_quadrant()関数により, 格子(x, y)を含む適切な子ノードを選択 */
    mark_mapsumpoint( tree->elem[find_quadrant( tree->center, x, y )],
		      map, x, y );
  } else {
    /* (x, y)に対応する格子のmapsumメンバを更新 */
    map->mapsum[tree->center.x/2][tree->center.y/2]++;
  }
}

/* 指定した格子に, スキャンデータのヒット回数を加算 */
void
set_maphitpoint( carmen_vascocore_map_t *map, int x, int y, float value  )
{
  mark_maphitpoint( &(map->qtree), map, 2*x, 2*y, value );
}

/* 指定した格子に, スキャンデータの参照回数を記録 */
void
set_mapsumpoint( carmen_vascocore_map_t *map, int x, int y )
{
  mark_mapsumpoint( &(map->qtree), map, 2*x, 2*y );
}

/* 使用されているノードの総数をカウント */
void
tree_list( carmen_vascocore_quad_tree_t *tree , int *ct )
{
  if ((tree->level)>0 ) {
    if (tree->elem[0]->inuse)
      tree_list( tree->elem[0], ct );
    if (tree->elem[1]->inuse)
      tree_list( tree->elem[1], ct );
    if (tree->elem[2]->inuse)
      tree_list( tree->elem[2], ct );
    if (tree->elem[3]->inuse)
      tree_list( tree->elem[3], ct );
  } else {
    (*ct)++;
  }
}

/* 4分木をクリア */
void
clear_local_treemap( carmen_vascocore_quad_tree_t *tree,
		     carmen_vascocore_map_t *map, int hk )
{
  int i,j;
  if ((tree->level)>0 ) {
    /* 再帰的にノードをクリア */
    if (tree->elem[0]->inuse)
      clear_local_treemap( tree->elem[0], map, hk );
    if (tree->elem[1]->inuse)
      clear_local_treemap( tree->elem[1], map, hk );
    if (tree->elem[2]->inuse)
      clear_local_treemap( tree->elem[2], map, hk );
    if (tree->elem[3]->inuse)
      clear_local_treemap( tree->elem[3], map, hk );
  } else {
    /* 葉ノードのクリア */
    /* 各ノードは, 領域の中心のインデックスを2倍した値を保持 */
    /* tree->center.x / 2 >= hk &&
     * tree->center.x / 2 < (map->mapsize.x - hk) &&
     * tree->center.y / 2 >= hk &&
     * tree->center.y / 2 < (map->mapsize.y - hk) */
    if ( (tree->center.x/2)>hk-1 && (tree->center.x/2)<map->mapsize.x-hk &&
	    (tree->center.y/2)>hk-1 && (tree->center.y/2)<map->mapsize.y-hk ) {
      /* (2 * hk + 1) * (2 * hk + 1)の領域をクリア */
      /* hkはガウスカーネルの半径 */
      for (i=(tree->center.x/2)-hk;i<=(tree->center.x/2)+hk;i++) {
        for (j=(tree->center.y/2)-hk;j<=(tree->center.y/2)+hk;j++) {
          map->maphit[i][j]  = 0;
          map->mapsum[i][j]  = 0.0;
          map->mapprob[i][j] = carmen_vascocore_settings.local_map_std_val;
          map->calc[i][j]    = carmen_vascocore_settings.local_map_std_val;
          map->updated[i][j] = UPDT_NOT;
        }
      }
    }
  }
  tree->inuse = FALSE;
}

/* スキャンデータから局所地図を作成 */
void
create_local_treemap( carmen_vascocore_map_t         * map,
		      carmen_vascocore_extd_laser_t    data,
		      carmen_move_t                    movement )
{
  int                         i;
  carmen_point_t              rpos;
  carmen_ivec2_t              start, end;
  carmen_vec2_t               lpos;

  /* 姿勢の相対変化を逆方向に変換 */
  rpos  = carmen_point_backwards_from_move( movement );
  /* 逆方向の変化量を, 格子のインデックスの差分に変換 */
  compute_map_pos_from_pos( rpos, *map, &start );

  for (i=0;i<data.numvalues;i++) {
    /* スキャン点と対応する座標を計算 */
    lpos = carmen_laser_point( rpos,  data.val[i], data.angle[i] );
    /* 座標を格子のインデックスに変換 */
    if (compute_map_pos_from_vec2( lpos, *map, &end )) {
      /* スキャン点が地図内に収まっていれば, 対応する格子の参照回数を更新 */
      if ( data.val[i]<carmen_vascocore_settings.local_map_max_range &&
          end.x>=0 && end.x<map->mapsize.x &&
          end.y>=0 && end.y<map->mapsize.y ) {
        set_mapsumpoint( map, end.x, end.y );
      }
    }
  }
}

/* スキャン点と対応する座標を計算 */
carmen_vec2_t
vascocore_compute_laser2d_coord( carmen_vascocore_extd_laser_t data, int i )
{
  return( carmen_laser_point( data.estpos, data.val[i], data.angle[i] ) );
}

/* 指定された距離だけ離れているか, または
 * 指定された角度だけ異なるかを判定 */
int
minimal_rpos_diff( carmen_point_t pos1, carmen_point_t pos2,
		   double pos_diff_min_dist,
		   double pos_diff_min_rot )
{
  carmen_vec2_t v1, v2;
  v1.x = pos1.x;    v1.y = pos1.y;
  v2.x = pos2.x;    v2.y = pos2.y;
  if ( carmen_vec_distance(v1,v2) > pos_diff_min_dist )
    return(TRUE);
  if ( carmen_orientation_diff(pos1.theta,pos2.theta) > pos_diff_min_rot )
    return(TRUE);
  return(FALSE);
}

/* 指定された距離だけ動いているか, または
 * 指定された角度だけ回転しているかを判定 */
int
minimal_rmove_diff( carmen_move_t move,
		    double pos_diff_min_dist,
		    double pos_diff_min_rot )
{
  carmen_vec2_t v1;
  v1.x = move.forward;    v1.y = move.sideward;
  if ( carmen_vec_length(v1) > pos_diff_min_dist )
    return(TRUE);
  if ( move.rotation > pos_diff_min_rot )
    return(TRUE);
  return(FALSE);
}

/* スキャンデータを囲むバウンディングボックスを計算 */
void
vascocore_compute_bbox( carmen_vascocore_extd_laser_t *data )
{
  int i;
  carmen_vec2_t min,max;
  min.x = DBL_MAX;     min.y = DBL_MAX;
  max.x = -DBL_MAX;    max.y = -DBL_MAX;
  /* スキャンデータの最小および最大の座標を計算 */
  for (i=0;i<data->numvalues;i++) {
    if (data->val[i]<carmen_vascocore_settings.bounding_box_max_range) {
      if (data->coord[i].x<min.x)
	      min.x = data->coord[i].x;
      if (data->coord[i].y<min.y)
	      min.y = data->coord[i].y;
      if (data->coord[i].x>max.x)
	      max.x = data->coord[i].x;
      if (data->coord[i].y>max.y)
	      max.y = data->coord[i].y;
    }
  }
  /* バウンディングボックスを指定された量だけ拡張 */
  min.x -= carmen_vascocore_settings.bounding_box_border;
  min.y -= carmen_vascocore_settings.bounding_box_border;
  max.x += carmen_vascocore_settings.bounding_box_border;
  max.y += carmen_vascocore_settings.bounding_box_border;
  /* バウンディングボックスを設定 */
  data->bbox.min = min;
  data->bbox.max = max;
}

/* 2つのバウンディングボックスの重なりを判定 */
int
intersect_bboxes( carmen_bbox_t box1, carmen_bbox_t box2 )
{
  /* max(box1.min.x, box2.min.x) < min(box1.max.x, box2.max.x) &&
   * max(box1.min.y, box2.min.y) < min(box1.max.y, box2.max.y) */
  if (box1.min.x<=box2.min.x) {
    /* box1.min.x is smaller that box2 */
    if (box1.max.x>box2.min.x) {
      /* intersection in x */
      if (box1.min.y<=box2.min.y) {
	/* box1.min.y is smaller that box2 */
	if (box1.max.y>box2.min.y) {
	  /* intersection in y */
	  return(1);
	} else {
	  return(0);
	}
      } else {
	/* box2.min.y is smaller that box1 */
	if (box2.max.y>=box1.min.y) {
	  /* intersection in y */
	  return(1);
	} else {
	  return(0);
	}
      }
    } else {
      return(0);
    }
  } else {
    /* box2.min.x is smaller that box1 */
    if (box2.max.x>=box1.min.x) {
      /* intersection in x */
      if (box1.min.y<=box2.min.y) {
	/* box1.min.y is smaller that box2 */
	if (box1.max.y>box2.min.y) {
	  /* intersection in y */
	  return(1);
	} else {
	  return(0);
	}
      } else {
	/* box2.min.y is smaller that box1 */
	if (box2.max.y>=box1.min.y) {
	  /* intersection in y */
	  return(1);
	} else {
	  return(0);
	}
      }
    } else {
      return(0);
    }
  }
  return(0);
}

