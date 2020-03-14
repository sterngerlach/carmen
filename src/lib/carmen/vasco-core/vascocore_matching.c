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

#include <sys/types.h>
#include <float.h>

#include "vasco.h"
#include "vascocore_intern.h"

/* 指定された格子の占有確率を取得 */
double
get_map_val( carmen_ivec2_t pos, carmen_vascocore_map_t map )
{
  if ( pos.x>=0 && pos.x<map.mapsize.x &&
       pos.y>=0 && pos.y<map.mapsize.y ) {
      return( EPSILON + (double) (map.mapprob[pos.x][pos.y]) );
  } else {
      return( EPSILON + carmen_vascocore_settings.local_map_std_val );
  }
  return(0.0);
}

/* 指定された格子の参照回数を取得 */
int
get_map_sum( carmen_ivec2_t pos, carmen_vascocore_map_t map )
{
  if ( pos.x>=0 && pos.x<map.mapsize.x &&
       pos.y>=0 && pos.y<map.mapsize.y ) {
      return( map.mapsum[pos.x][pos.y] );
  } else {
    return( 0 );
  }
  return(0);
}

/* 指定された格子にレーザスキャンが当たった回数を取得 */
float
get_map_hit( carmen_ivec2_t pos, carmen_vascocore_map_t map )
{
  if ( pos.x>=0 && pos.x<map.mapsize.x &&
       pos.y>=0 && pos.y<map.mapsize.y ) {
      return( map.maphit[pos.x][pos.y] );
  } else {
    return( 0 );
  }
  return(0);
}


#define LOG_EPSILON                -100

/* ある移動量の下で, 指定された移動量が得られる確率(対数)を計算 */
double
probability_between_moves_old( carmen_move_t move1, carmen_move_t move2 )
{
  double val1, val2, val3;

  val1   = fabs( move2.forward - move1.forward );
  val2   = fabs( move2.sideward - move1.sideward );
  val3   = fabs( carmen_orientation_diff( move1.rotation,
					  move2.rotation ) );

  if (val1>carmen_vascocore_settings.motion_model_forward) {
    return( LOG_EPSILON );
  }
  if (val2>carmen_vascocore_settings.motion_model_sideward) {
    return( LOG_EPSILON );
  }
  if (val3>carmen_vascocore_settings.motion_model_rotation) {
    return( LOG_EPSILON );
  }

  return(  log( (carmen_vascocore_settings.motion_model_forward - val1)/
		carmen_vascocore_settings.motion_model_forward ) +
	   log( (carmen_vascocore_settings.motion_model_sideward - val2)/
		carmen_vascocore_settings.motion_model_sideward ) +
	   log( (carmen_vascocore_settings.motion_model_rotation - val3)/
		carmen_vascocore_settings.motion_model_rotation ) );
}

/* ある移動量の下で, 指定された移動量が得られる確率(対数)を計算 */
double
probability_between_moves( carmen_move_t move1, carmen_move_t move2 )
{
  double sum = 0.0;
  sum += log( EPSILON +
	      carmen_gauss( fabs( move1.forward-move2.forward ),
			    0, carmen_vascocore_settings.motion_model_forward ));
  sum += log( EPSILON +
	      carmen_gauss( fabs( move1.sideward-move2.sideward ),
			    0, carmen_vascocore_settings.motion_model_sideward));
  sum += log( EPSILON +
	      carmen_gauss( fabs( carmen_orientation_diff( move1.rotation,
							   move2.rotation )),
			    0, carmen_vascocore_settings.motion_model_rotation ));
  return( sum );
}

/* 予測値と観測値との差分から誤差値を計算 */
double
error( double val, double expect )
{
  double sigma1=40.0;
  double sigma2=80.0;

  /* valがexpectよりもsigma2以上離れている場合は飽和 */
  /* valがexpect以上かつ, valがexpectよりもsigma1以上離れている場合
   * (valがexpectよりもsigma1以上大きい場合)は飽和 */
  /* val <= expect - sigma2 または val >= expect + sigma2 であれば 1.0
   * val >= expect + sigma1 であれば 1.0
   * expect - sigma2 < val <= 0 であれば (val - expect) / sigma2 (1.0から0.0)
   * 0 < val <= expect + sigma1 であれば (expect - val) / sigma2 (0.0から0.5) */
  if ( fabs(val-expect)<sigma2 &&
       (val<expect || fabs(val-expect)<sigma1) ) {
    /* return fabs(val - expect) / sigma2; */
    return( 1.0- (sigma2-fabs(val-expect))/sigma2);
  } else {
    /* valが離れすぎている場合は誤差値を飽和 */
    return( 1.0 );
  }
}

/* 2つの整数ベクトル間の距離を計算 */
double
ivector2_distance( carmen_ivec2_t p1, carmen_ivec2_t p2 )
{
  return sqrt( (p1.x-p2.x)*(p1.x-p2.x) +
	       (p1.y-p2.y)*(p1.y-p2.y) );
}

/* 測定距離の予測値と観測値から, 観測値が得られる確率を計算 */
double
compute_beam_log_prob( double expected, double measured )
{
  double val, d = fabs(expected-measured); /* dist in cm */

  /* measuredが最大値に近い場合 */
  if (measured>0.95*carmen_vascocore_settings.local_map_max_range)
    return(log(0.01));
  if (d>carmen_vascocore_settings.local_map_max_range)
    d = carmen_vascocore_settings.local_map_max_range;
  if (d<200.0) {
    val = (220.0-d)/220.0;
  } else {
    /* expectedとmeasuredの差分が大き過ぎる場合 */
    val = 0.01;
  }
  return(log(val));

}

/* 動作モデルで重み付けされた計測尤度を計算 */
double
probability_with_move( carmen_vascocore_map_t           map,
		       carmen_vascocore_extd_laser_t    data,
		       carmen_move_t                    move,
		       carmen_move_t                    odo_move,
		       double                         * laserprob )
{
  int                   i;
  double                bprob, prob = 0.0;
  carmen_vec2_t         pt;
  carmen_ivec2_t        mvec;
  carmen_point_t        rpos;

  rpos = carmen_point_from_move( move );
  /* 各スキャンデータに対する処理 */
  for (i=0;i<data.numvalues;i++) {
    /* 距離値が妥当である(最大値未満である)場合 */
    if (data.val[i]<carmen_vascocore_settings.local_map_max_range) {
      /* 距離値と角度から局所地図上での座標を計算 */
      pt = carmen_laser_point( rpos, data.val[i], data.angle[i] );
      /* 局所地図上での座標を格子のインデックスに変換 */
      compute_map_pos_from_vec2( pt, map, &mvec );
      /* 格子の占有確率の対数値を算出 */
      bprob = log(get_map_val( mvec, map ));
    } else {
      /* 距離値が無効であれば, デフォルトの尤度値に設定 */
      bprob = log(carmen_vascocore_settings.local_map_std_val);
    }
    /* 対数尤度に, 各スキャンデータに対する尤度値を加算 */
    prob += bprob;
  }

  *laserprob = prob;
  /* 動作モデルの対数値を計算し尤度に加算 */
  if (carmen_vascocore_settings.local_map_use_odometry)
    prob += probability_between_moves( move, odo_move );

  return( prob );
}

/* 動作モデルで重み付けされた計測尤度を計算 */
double
probability_with_move_new( carmen_vascocore_map_t           map,
			   carmen_vascocore_extd_laser_t    data,
			   carmen_move_t                    move,
			   carmen_move_t                    odo_move,
			   double                         * laserprob )
{
  int                   i;
  double                bprob, prob = 0.0;
  carmen_vec2_t         pt;
  carmen_ivec2_t        mvec;
  carmen_point_t        rpos;

  rpos = carmen_point_from_move( move );

  for (i=0;i<data.numvalues;i++) {
    if (data.val[i]<carmen_vascocore_settings.local_map_max_range) {
      pt = carmen_laser_point( rpos, data.val[i], data.angle[i] );
      compute_map_pos_from_vec2( pt, map, &mvec );
      bprob = log(get_map_val( mvec, map ));
    } else {
      bprob = log(carmen_vascocore_settings.local_map_std_val);
    }
    prob += bprob;
  }

  *laserprob = prob;
  if (carmen_vascocore_settings.local_map_use_odometry)
    prob += probability_between_moves( move, odo_move );

  return( prob );
}

/* 指定されたセンサ位置の下で, 動作モデルで重み付けされた計測尤度を計算 */
double
probability_with_pos( carmen_vascocore_map_t           map,
		      carmen_vascocore_extd_laser_t    data,
		      carmen_point_t                   rpos,
		      carmen_move_t                    move,
		      carmen_move_t                    odo_move )
{
  int              i;
  carmen_vec2_t          pt;
  carmen_ivec2_t         mvec;
  double           bprob, prob = 0.0;

  for (i=0;i<data.numvalues;i++) {
    if (data.val[i]<carmen_vascocore_settings.local_map_max_range) {
      pt = carmen_laser_point( rpos, data.val[i], data.angle[i] );
      compute_map_pos_from_vec2( pt, map, &mvec );
      bprob = log(get_map_val( mvec, map ));
    } else {
      bprob = log(carmen_vascocore_settings.local_map_std_val);
    }
    prob += bprob;
  }
  prob += probability_between_moves( move, odo_move );

  return( prob );
}

/* 指定された方向に少しだけ移動した移動量を計算 */
carmen_move_t
compute_test_move( carmen_move_t smove, int nummove, int stepsize )
{
  carmen_move_t move = smove;
  double div  = pow( 2, stepsize);
  switch( nummove ) {
  case 0:
    move.rotation += ( carmen_vascocore_settings.pos_corr_step_size_rotation / div );
    break;
  case 1:
    move.rotation -= ( carmen_vascocore_settings.pos_corr_step_size_rotation / div );
    break;
  case 2:
    move.sideward += ( carmen_vascocore_settings.pos_corr_step_size_sideward / div );
    break;
  case 3:
    move.sideward -= ( carmen_vascocore_settings.pos_corr_step_size_sideward / div );
    break;
  case 4:
    move.forward  += ( carmen_vascocore_settings.pos_corr_step_size_forward / div) ;
    break;
  case 5:
    move.forward  -= ( carmen_vascocore_settings.pos_corr_step_size_forward / div );
    break;
  default:
    break;
  }
  return(move);
}

/* 計測尤度を最大にするような移動量を山登り法により計算
 * 地図とスキャンデータとの重ね合わせを実行 */
carmen_move_t
fit_data_in_local_map( carmen_vascocore_map_t            map,
		       carmen_vascocore_extd_laser_t     data,
		       carmen_move_t                     movement )
{
  int i, l;
  int fitting = TRUE;
  carmen_move_t bmove, pmove, tmove;
  double bprob, pprob, prob, laserprob;

  int loop = 0, adjusting = TRUE;

  /* bprob = probability_with_move( map, data, movement, movement, &laserprob ); */

  /* 計測尤度と姿勢を初期化 */
  pmove = bmove = movement;
  bprob = -FLT_MAX;

  l = 0;

  while( adjusting ) {

    loop    = 0;
    fitting = TRUE;

    while( fitting ) {

      pprob = bprob;
      /* 前後左右に移動, あるいは左右に回転したときの, 計測尤度を計算 */
      for (i=0; i<6; i++) {
        tmove = compute_test_move( bmove, i, loop );
        prob = probability_with_move( map, data, tmove,
                    movement, &laserprob  );
        /* 最も計測尤度が高くなるときの姿勢をpmove, 計測尤度をpprobに保持 */
        if ( prob>pprob) {
          pmove  = tmove;
          pprob = prob;
        }
      }

      /* 移動によって計測尤度が改善した場合は, 計測尤度の最善値を更新 */
      if (pprob-bprob>EPSILON) {
        bmove  = pmove;
        bprob = pprob;
      } else if (loop<carmen_vascocore_settings.pos_corr_step_size_loop) {
        /* 計測尤度が改善しなくなった場合はループ回数を加算
         * compute_test_moveの第3引数が変化するので, 姿勢の変化量が減衰する */
	      loop++;
      } else {
        /* ループ回数が指定された回数を超えたら終了 */
	      fitting = FALSE;
      }
    } /* end while( fitting ) */

    l++;

    if (l>0)
      adjusting = FALSE;

  } /* end while( adjust ) */

  return(bmove);
}
