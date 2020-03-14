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

#ifndef VASCOCORE_H
#define VASCOCORE_H

#define EPSILON                          0.00000001
#define MAX_NUM_LASER_VALUES           401

#define UPDT_NOT                         0
#define UPDT_X                           1
#define UPDT_Y                           2

typedef struct carmen_vascocore_quad_tree_t {
  
  /* 右上, 左上, 左下, 右下の子ノード */
  struct carmen_vascocore_quad_tree_t  * elem[4];
  /* 領域の中心のインデックスを, 2倍にした値 */
  carmen_svec2_t                         center;
  /* ノードの深さ */
  unsigned char                          level;
  /* ノードの使用状態(0: 未使用, 1: 使用) */
  unsigned char                          inuse;

} carmen_vascocore_quad_tree_t;

/* 占有格子地図を表す構造体 */
typedef struct {
  
  /* 地図をカバーする4分木 */
  carmen_vascocore_quad_tree_t           qtree;
  carmen_point_t                         offset;
  /* 各格子のサイズ */
  double                                 resolution;
  /* 畳み込みの中間結果(calc)の状態 */
  unsigned char                       ** updated;
  /* 各格子にスキャンデータがヒットした(当たった)回数 */
  float                               ** maphit;
  /* 各格子にスキャンデータがヒット(当たった)および
   * ミスした(格子上をスキャンデータが通過した)回数の合計
   * 各格子がスキャンによって参照された回数 */
  short                               ** mapsum;
  /* 各格子の占有確率(ヒットした回数と参照回数との比率 */
  float                               ** mapprob;
  /* 畳み込みの中間結果の値 */
  float                               ** calc;
  /* 地図の半径(横方向と縦方向の格子数の半分) */
  carmen_ivec2_t                         mapsize;
  /* 地図の中心の座標 */
  carmen_vec2_t                          center;

} carmen_vascocore_map_t, *carmen_vascocore_map_p;

/* スキャンデータを表す構造体 */
typedef struct {

  double                                 time;
  /* スキャンデータ取得時のセンサの位置 */
  carmen_point_t                         estpos;
  double                                 fov;
  int                                    numvalues;
  /* センサの中心からの距離 */
  double                               * val;
  /* センサの中心からの角度 */
  double                               * angle;
  /* 直交座標系での表現 */
  carmen_vec2_t                        * coord;
  /* スキャンデータを囲むバウンディングボックス */
  carmen_bbox_t                          bbox;

} carmen_vascocore_extd_laser_t;

/* スキャンデータの履歴を表す構造体 */
typedef struct {
  
  int                                    started;
  int                                    length;
  int                                    ptr;
  carmen_vascocore_extd_laser_t        * data;
  
} carmen_vascocore_history_t;

extern carmen_vascocore_param_t      carmen_vascocore_settings;
extern carmen_vascocore_map_t        carmen_vascocore_map;
extern carmen_vascocore_history_t    carmen_vascocore_history;

void   initialize_map( carmen_vascocore_map_t * map,
		       int sx, int sy, int center_x, int center_y,
		       double resolution, carmen_point_t start );

void   clear_local_treemap( carmen_vascocore_quad_tree_t *tree,
			    carmen_vascocore_map_t *map, int hk );

void   create_local_treemap( carmen_vascocore_map_t * map,
			     carmen_vascocore_extd_laser_t data,
			     carmen_move_t movement );

void   convolve_treemap( carmen_vascocore_map_t *map );

int    intersect_bboxes( carmen_bbox_t box1, carmen_bbox_t box2 );

carmen_move_t
fit_data_in_local_map( carmen_vascocore_map_t map,
		       carmen_vascocore_extd_laser_t data,
		       carmen_move_t movement );

void   vascocore_compute_bbox( carmen_vascocore_extd_laser_t *data );

carmen_vec2_t
vascocore_compute_laser2d_coord( carmen_vascocore_extd_laser_t data, int i );

int     compute_map_pos_from_vec2( carmen_vec2_t pos,
				   carmen_vascocore_map_t map,
				   carmen_ivec2_t *v );

void           vascocore_graphics_map( void );

void           vascocore_graphics_update( void );

void           vascocore_graphics_init( int argc, char *argv[] );

void           vascocore_graphics_scan( carmen_vascocore_extd_laser_t scan );

#endif
