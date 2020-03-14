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

#include "vasco.h"
#include "vascocore_intern.h"

/* 局所地図をクリア */
void
clear_local_map( void )
{
  clear_local_treemap( &(carmen_vascocore_map.qtree), &carmen_vascocore_map,
		       (carmen_vascocore_settings.local_map_kernel_len-1)/2 );
}

/* スキャンデータから局所地図を作成 */
void
create_local_map( carmen_vascocore_extd_laser_t data, carmen_move_t move )
{
  create_local_treemap( &carmen_vascocore_map, data, move );
}

/* 局所地図を畳み込み */
void
convolve_map( void )
{
  convolve_treemap( &carmen_vascocore_map );
}

/* 計測尤度を最大化するときの移動量を計算
 * 地図とスキャンデータとの重ね合わせを実行 */
carmen_move_t
find_best_move( carmen_vascocore_extd_laser_t data, carmen_move_t move )
{
  return( fit_data_in_local_map( carmen_vascocore_map, data, move ) );
}

/* スキャンデータのメッセージをVascoモジュールの構造体にコピー */
void
vascocore_copy_scan( carmen_laser_laser_message     scan,
		     carmen_vascocore_extd_laser_t *data )
{
  int       i;
  double    delta,start;

  data->time = scan.timestamp;

  data->fov =  scan.config.fov;
  delta = scan.config.angular_resolution;
  start = scan.config.start_angle;

  for (i=0;i<scan.num_readings; i++) {
    data->val[i] = (double) scan.range[i];
    data->angle[i] = start + i*delta;
  }
  data->numvalues = scan.num_readings;
}

/* スキャンデータを格納する履歴のインデックスを決定 */
int
hpos( int pos )
{
  return(pos % carmen_vascocore_settings.local_map_history_length);
}

/* スキャンマッチングを実行 */
carmen_point_t
vascocore_scan_match( carmen_laser_laser_message scan, carmen_point_t pos )
{
  /* オドメトリにより得られたロボットの姿勢 */
  static carmen_point_t  lastpos;
  carmen_point_t         estpos, centerpos, *histpos, nopos = {0.0, 0.0, 0.0};
  carmen_move_t          estmove, bestmove, move, nullmove = {0.0, 0.0, 0.0};
  int                    i, h, hp, hps, p, ctr=0;

  /* スキャンデータを格納する履歴のインデックスを計算 */
  p = hpos( carmen_vascocore_history.ptr );

  /* 現在のスキャンデータを履歴に保持 */
  vascocore_copy_scan( scan, &(carmen_vascocore_history.data[p]) );

  /* 最初のスキャンデータであれば, 初期姿勢は原点とする */
  if (!carmen_vascocore_history.started) {
    /* THE FIRST SCAN WILL BE MAPPED TO 0/0 */
    estpos.x      = 0.0;
    estpos.y      = 0.0;
    estpos.theta  = 0.0;
    if (carmen_vascocore_settings.verbose) {
      fprintf( stderr, "***************************************\n" );
      fprintf( stderr, "* first scan ...\n" );
    }
    carmen_vascocore_history.started = TRUE;
  } else {
    /* startedがTRUEであるので, 少なくとも1つのスキャンデータは履歴に保持されている */
    /* 直前の姿勢からの変化(オドメトリ)を計算 */
    estmove = carmen_move_between_points( lastpos, pos );
    /* 局所地図のクリア */
    clear_local_map();
    /* 直前のスキャンデータの格納先のインデックスを取得 */
    hps = hpos(carmen_vascocore_history.ptr-1);
    /* 直前のスキャンデータ取得時の姿勢を取得 */
    /* 局所地図の原点と対応するロボットの座標 */
    centerpos = carmen_vascocore_history.data[hps].estpos;
    /* 直前のスキャンデータから局所地図を作成 */
    create_local_map( carmen_vascocore_history.data[hps], nullmove );
    /* 直前のロボットの姿勢を取得 */
    histpos = &carmen_vascocore_history.data[hps].estpos;
    /* 指定されたステップ数(local_map_history_length)だけ遡る */
    /* 最大local_map_max_used_history個のスキャンデータを元に部分地図を作成 */
    for (  h=carmen_vascocore_history.ptr-2;
	   h>=0 &&
	   h>(carmen_vascocore_history.ptr-
	      carmen_vascocore_settings.local_map_history_length) &&
	   ctr<carmen_vascocore_settings.local_map_max_used_history;
	   h--) {
      /* 更に前に遡ったときのスキャンデータの履歴のインデックス */
      hp = hpos(h);
      /* 直前のスキャンデータとバウンディングボックスが重なり合っている場合 */
      /* 最低local_map_use_last_scans個分のスキャンデータから局所地図を作成 */
      /* 現在みているスキャンデータと, その1つ手前のスキャンデータが取得されたときの,
       * 2つの姿勢間の距離が一定以上であれば利用 */
      if ( intersect_bboxes( carmen_vascocore_history.data[hps].bbox,
			     carmen_vascocore_history.data[hp].bbox ) &&
          ( ( (carmen_vascocore_history.ptr-1-h) <
              carmen_vascocore_settings.local_map_use_last_scans ) ||
            (carmen_point_dist( *histpos,
              carmen_vascocore_history.data[hp].estpos ) >
              carmen_vascocore_settings.local_map_min_bbox_distance)) ) {
        /* 直前のスキャンデータが取得されたときの姿勢と,
         * 更に遡ったスキャンデータが取得されたときの姿勢との相対変化を計算 */
        move = carmen_move_between_points( carmen_vascocore_history.data[hp].estpos,
                  centerpos );
        /* 姿勢の相対変化を考慮し, 遡ったスキャンデータを基に局所地図を更新 */
        create_local_map( carmen_vascocore_history.data[hp], move );
        /* 遡ったときの姿勢を保持 */
        histpos = &carmen_vascocore_history.data[hp].estpos;
        ctr++;
      }
    }
    /* COMPUTE AND CONVOLVE LOCAL MAP */
    /* 局所地図を畳み込み */
    convolve_map();
    if (carmen_vascocore_settings.verbose) {
      fprintf( stderr, "***************************************\n" );
      fprintf( stderr, "using %d scans\n", ctr );
      fprintf( stderr, "***************************************\n" );
      fprintf( stderr, "estimated movement    %.4f %.4f %.4f\n",
	      estmove.forward, estmove.sideward, estmove.rotation );
    }
    /* 局所地図とスキャンデータとの重ね合わせにより移動量を改良 */
    bestmove = find_best_move( carmen_vascocore_history.data[p], estmove );
    if (carmen_vascocore_settings.verbose) {
      fprintf( stderr, "best movement         %.4f %.4f %.4f\n",
	      bestmove.forward, bestmove.sideward, bestmove.rotation );
    }
    /* 改良された移動量を, 直前のスキャンデータが取得されたときの姿勢に反映 */
    /* 現在のスキャンデータの取得時の姿勢を計算 */
    estpos = carmen_point_with_move( centerpos, bestmove );
  }

  /* SAVE THE COMPUTED POSITION AND COORDS */
  /* ここでestposを設定するのは, vascocore_compute_laser2d_coord関数で使用するため */
  carmen_vascocore_history.data[p].estpos =
    carmen_point_with_move( nopos, bestmove );
  for (i=0; i<carmen_vascocore_history.data[p].numvalues; i++) {
    /* 各スキャンデータの距離と角度から直交座標を計算 */
    /* スキャンマッチングでは結局使われない */
    carmen_vascocore_history.data[p].coord[i] =
      vascocore_compute_laser2d_coord( carmen_vascocore_history.data[p], i );
  }
  /* ここで最終的なestposを設定 */
  carmen_vascocore_history.data[p].estpos = estpos;
  /* 現在のスキャンデータのバウンディングボックスを計算 */
  vascocore_compute_bbox( &carmen_vascocore_history.data[p] );

  /* オドメトリによる姿勢を更新 */
  lastpos = pos;
  /* 次のスキャンデータを格納すべき, 履歴のインデックスを更新 */
  carmen_vascocore_history.ptr++;

  return(estpos);

}

/* スキャンデータをVascoモジュールの構造体にコピー */
void
vascocore_copy_scan_general(int num_readings, float *range, float *angle,
			    double fov,
			    carmen_vascocore_extd_laser_t *data)
{
  int       i;

  data->time = 0;
  data->numvalues = num_readings;
  for(i = 0; i < num_readings; i++) {
    data->val[i] = range[i];
    data->angle[i] = angle[i];
  }
  data->fov = fov;
}

/* スキャンマッチングを実行 */
carmen_point_t
vascocore_scan_match_general(int num_readings, float *range, float *angle,
			     double fov,
			     carmen_point_t pos, int first)
{
  static carmen_point_t  lastpos;
  carmen_point_t         estpos, centerpos, *histpos, nopos = {0.0, 0.0, 0.0};
  carmen_move_t          estmove, bestmove, move, nullmove = {0.0, 0.0, 0.0};
  int                    i, h, hp, hps, p, ctr=0;

  p = hpos( carmen_vascocore_history.ptr );

  vascocore_copy_scan_general(num_readings, range,
			      angle,
			      fov,
			      &(carmen_vascocore_history.data[p]));

  if(first) {
    /* THE FIRST SCAN WILL BE MAPPED TO 0/0 */
    estpos.x      = 0.0;
    estpos.y      = 0.0;
    estpos.theta  = 0.0;
    if (carmen_vascocore_settings.verbose) {
      fprintf( stderr, "***************************************\n" );
      fprintf( stderr, "* first scan ...\n" );
    }
  } else {
    estmove = carmen_move_between_points( lastpos, pos );
    clear_local_map();
    hps = hpos(carmen_vascocore_history.ptr-1);
    centerpos = carmen_vascocore_history.data[hps].estpos;
    /* CREATE LOCAL MAP FROM HISTORY */
    create_local_map( carmen_vascocore_history.data[hps], nullmove );
    histpos = &carmen_vascocore_history.data[hps].estpos;
    for (  h=carmen_vascocore_history.ptr-2;
	   h>=0 &&
	   h>(carmen_vascocore_history.ptr-
	      carmen_vascocore_settings.local_map_history_length) &&
	   ctr<carmen_vascocore_settings.local_map_max_used_history;
	   h--) {
      hp = hpos(h);
      if ( intersect_bboxes( carmen_vascocore_history.data[hps].bbox,
			     carmen_vascocore_history.data[hp].bbox ) &&
            carmen_point_dist( *histpos,
                    carmen_vascocore_history.data[hp].estpos ) >
            carmen_vascocore_settings.local_map_min_bbox_distance ) {
        move = carmen_move_between_points( carmen_vascocore_history.data[hp].estpos,
                  centerpos );
        create_local_map( carmen_vascocore_history.data[hp], move );
        histpos = &carmen_vascocore_history.data[hp].estpos;
        ctr++;
      }
    }
    /* COMPUTE AND CONVOLVE LOCAL MAP */
    convolve_map();
    if (carmen_vascocore_settings.verbose) {
      fprintf( stderr, "***************************************\n" );
      fprintf( stderr, "using %d scans\n", ctr );
      fprintf( stderr, "***************************************\n" );
      fprintf( stderr, "estimated movment    %.4f %.4f %.4f\n",
	      estmove.forward, estmove.sideward, estmove.rotation );
    }
    bestmove = find_best_move( carmen_vascocore_history.data[p], estmove );
    if (carmen_vascocore_settings.verbose) {
      fprintf( stderr, "best movment         %.4f %.4f %.4f\n",
	      bestmove.forward, bestmove.sideward, bestmove.rotation );
    }
    estpos = carmen_point_with_move( centerpos, bestmove );
  }

  /* SAVE THE COMPUTED POSITION AND COORDS */
  carmen_vascocore_history.data[p].estpos =
    carmen_point_with_move( nopos, bestmove );
  for (i=0; i<carmen_vascocore_history.data[p].numvalues; i++) {
    carmen_vascocore_history.data[p].coord[i] =
      vascocore_compute_laser2d_coord( carmen_vascocore_history.data[p], i );
  }
  carmen_vascocore_history.data[p].estpos = estpos;
  vascocore_compute_bbox( &carmen_vascocore_history.data[p] );

  lastpos = pos;
  carmen_vascocore_history.ptr++;

  return(estpos);

}

