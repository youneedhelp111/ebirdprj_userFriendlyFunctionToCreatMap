## ============================================================================
## map_local() —— 用本地 eBird 基本数据集 (EBD/SED) 从零建模, 一键出相对丰度图
## ----------------------------------------------------------------------------
## 你只需要在工程 /data 目录放好 eBird Custom Download 的两个纯文本:
##   * 清单文件 SED:  文件名里带 "sampling",  例如 ebd_..._sampling_relJul-2026.txt
##   * 观测文件 EBD:  文件名以 "ebd" 开头但不带 "sampling", 例如 ebd_..._relJul-2026.txt
##   (在 eBird Custom Download 里勾选 "Include sampling event data" 才会同时给这两个)
## 函数自动完成:
##   读取/过滤 -> zero-fill(检测/非检测) -> 努力量派生变量 -> 时空子采样
##   -> 在线下载环境栅格(土地覆盖+高程, 免密钥) -> 提取环境变量
##   -> hurdle 双随机森林(遭遇率 RF + SCAM 校准 + 计数 RF)
##   -> 预测网格 -> 相对丰度 = 校准遭遇率 × 计数 -> 出图
##
## 依赖: dplyr, sf, terra, lubridate; 首次会自动从 CRAN 安装
##       ranger, scam, mccf1, elevatr, exactextractr, readr, hms, tidyr,
##       jsonlite, httr, rnaturalearth, rnaturalearthdata, fields
## ============================================================================


## ============================================================================

## ----------------------------------------------------------------------------
## 0. 依赖检查 / 安装
## ----------------------------------------------------------------------------
.check_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE),
                          logical(1))]
  if (length(missing)) {
    stop("缺少以下 R 包, 请先安装: ",
         paste(missing, collapse = ", "),
         "\n安装命令: install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), "))",
         call. = FALSE)
  }
}

## custom 路径专用: 缺包就自动从 CRAN 安装
.ensure_pkgs <- function(pkgs, verbose = TRUE) {
  missing <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE),
                          logical(1))]
  if (length(missing)) {
    if (verbose) message(">>> 首次运行 custom, 自动安装依赖包: ",
                         paste(missing, collapse = ", "))
    utils::install.packages(missing, repos = "https://cloud.r-project.org",
                            quiet = TRUE)
  }
  still <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE),
                        logical(1))]
  if (length(still)) stop("这些包安装失败, 请手动 install.packages: ",
                          paste(still, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

## NULL 值兜底运算符: a 为 NULL/长度0 时返回 b
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## ============================================================================
## 1. 主函数
## ============================================================================
#
##' 一键生成 eBird 相对丰度地图
##'
##' @param species         物种英文名(common name)或学名; 如 "Wood Thrush"。
##'                        (与 EBD 里 common_name 一致), 如 "Wood Thrush"。
##' @param season          时段: breeding / nonbreeding / prebreeding_migration /
##'                        postbreeding_migration / year_round / weekly / custom。
##' @param region          研究区。支持:
##'                        list(type="state",  name="Georgia", country_iso="US")
##'                        list(type="country",name="United States of America")
##'                        list(type="conus")                       (美国本土48州)
##'                        list(type="none")                        (用数据全范围)
##'                        或直接给一个 sf 多边形对象。
##' @param week            season="weekly" 时, 第几周 (1-52)。
##' @param start_date,end_date  season="custom" 时的起止日期。
##' @param data_dir        [custom] EBD/SED 所在目录, 默认 "data"。
##' @param ebd_file,sed_file    [custom] 可显式指定观测/清单文件名; NULL 自动探测。
##' @param year_range      [custom] 建模年份范围(整数向量); NULL=自动用数据全部年份。
##' @param grid_res        [custom] 预测网格边长(米), 默认 3000。
##' @param elev_z          [custom] elevatr 高程缩放级别 (7≈1km, 8≈600m)。
##' @param lc_agg          [custom] 土地覆盖聚合分辨率(米), 默认 90, 越大越快。
##' @param sample_radius   [custom] 环境变量提取的圆形邻域半径(米), 默认 1500(直径3km)。
##' @param compute_ed      [custom] 是否计算土地覆盖边缘密度(很慢, 默认 FALSE)。
##' @param test_fraction   [custom] 随机留出测试集比例, 默认 0.2。
##' @param ss_cellsize     [custom] 时空子采样空间网格边长(米), 默认 3000。
##' @param ss_weeks        [custom] 时空子采样时间窗(周), 默认 1。
##' @param hours_peak      [custom] 预测用的最佳观测开始时刻(小时); NULL=自动估计。
##' @param n_trees         [custom] 随机森林树数, 默认 500。
##' @param env_cache_dir   [custom] 环境栅格缓存目录, 默认 "data-raw/env-cache"。
##' @param seed            [custom] 随机种子, 默认 1。
##' @param out_dir         输出目录, 默认 "output"。
##' @param out_prefix      文件名前缀; NULL 自动生成。
##' @param map_title       地图标题; NULL 自动生成。
##' @param crs             投影; "auto"=按区域质心选 LAEA。
##' @param verbose         是否打印进度, 默认 TRUE。
##' @param do_plot         是否立即出图, 默认 TRUE。
##' @param ...             其它预留参数。
##'
##' @return 一个 S3 对象 "rel_abundance":
##'   r       SpatRaster (custom 为四层: in_range/encounter_rate/count/abundance,
##'                        出图与 print 默认用其中的 abundance 层)
##'   region  sf 多边形; species, season, call ...
##' @export
map_local <- function(species,
                     season = c("breeding", "nonbreeding",
                                "prebreeding_migration",
                                "postbreeding_migration",
                                "year_round",
                                "weekly", "custom"),
                     region = list(type = "state", name = "Georgia",
                                   country_iso = "US"),
                     week = NULL,
                     start_date = NULL, end_date = NULL,
                     ## ---- 本地建模专用参数 ----
                     data_dir       = "data",
                     ebd_file       = NULL,
                     sed_file       = NULL,
                     year_range             = NULL,
                     grid_res       = 3000,
                     elev_z         = 7,
                     lc_agg         = 90,
                     sample_radius  = 1500,
                     compute_ed     = FALSE,
                     test_fraction  = 0.2,
                     ss_cellsize    = 3000,
                     ss_weeks       = 1,
                     hours_peak     = NULL,
                     n_trees        = 500,
                     env_cache_dir  = "data-raw/env-cache",
                     seed           = 1,
                     smooth         = TRUE,
                     ## ---- 通用参数 ----
                     out_dir = "output",
                     out_prefix = NULL,
                     map_title = NULL,
                     crs = "auto",
                     verbose = TRUE,
                     do_plot = TRUE,
                     ...) {

  season <- match.arg(season)

  ## 基础依赖: 检查通过后 attach, 保证 source() 后即使没手动 library 也能直接调用
  .check_pkgs(c("dplyr", "sf", "terra", "lubridate"))
  suppressPackageStartupMessages({
    library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
    library(sf, quietly = TRUE, warn.conflicts = FALSE)
    library(terra, quietly = TRUE, warn.conflicts = FALSE)
    library(lubridate, quietly = TRUE, warn.conflicts = FALSE)
  })

  ## 解析研究区为 sf 多边形 (WGS84)
  region_sf <- .resolve_region_ne(region)

  ## 自动地图标题
  if (is.null(map_title)) {
    map_title <- paste0(species, " — ", .season_label(season, start_date,
                                                      end_date, week))
  }

  res <- .run_custom(species = species, season = season,
                     region_sf = region_sf,
                     week = week, start_date = start_date,
                     end_date = end_date,
                     data_dir = data_dir, ebd_file = ebd_file,
                     sed_file = sed_file, year_range = year_range,
                     grid_res = grid_res, elev_z = elev_z,
                     lc_agg = lc_agg, sample_radius = sample_radius,
                     compute_ed = compute_ed,
                     test_fraction = test_fraction,
                     ss_cellsize = ss_cellsize, ss_weeks = ss_weeks,
                     hours_peak = hours_peak, n_trees = n_trees,
                     env_cache_dir = env_cache_dir, seed = seed,
                     smooth = smooth,
                     crs = crs, out_dir = out_dir, out_prefix = out_prefix,
                     map_title = map_title, do_plot = do_plot,
                     verbose = verbose, ...)
  invisible(res)
}

## ============================================================================
## 2. 从 EBD/SED 清单数据从零建模 (hurdle 模型, 复刻官方教程)
## ============================================================================

## ============================================================================
## 3. 路径 B: 从 EBD/SED 清单数据从零建模 (hurdle 模型, 复刻官方教程)
## ============================================================================
.run_custom <- function(species, season, region_sf, week, start_date, end_date,
                        data_dir, ebd_file, sed_file, year_range,
                        grid_res, elev_z, lc_agg, sample_radius, compute_ed,
                        test_fraction, ss_cellsize, ss_weeks,
                        hours_peak, n_trees, env_cache_dir, seed,
                        smooth, crs, out_dir, out_prefix, map_title,
                        do_plot, verbose, ...) {

  .ensure_pkgs(c("ranger", "scam", "mccf1", "elevatr", "exactextractr",
                 "readr", "hms", "tidyr", "jsonlite", "httr",
                 "rnaturalearth", "rnaturalearthdata", "fields"),
               verbose = verbose)
  suppressPackageStartupMessages({
    library(dplyr); library(sf); library(terra); library(lubridate)
    library(readr); library(tidyr); library(ranger); library(scam)
    library(mccf1); library(exactextractr)
  })
  set.seed(seed)

  if (verbose) message(">>> [1/8] 解析时段、定位数据 ...")
  yr_place <- if (is.null(year_range)) 2023 else max(year_range)
  sinfo <- .custom_season_info(season, week, start_date, end_date, yr_place)
  season_used <- season
  auto_region <- isTRUE(attr(region_sf, "auto_region"))

  global_bb <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -85, ymax = 85),
                           crs = 4326)
  read_bb <- if (auto_region) global_bb
             else sf::st_bbox(sf::st_transform(region_sf, 4326)) +
                    c(-.1, -.1, .1, .1)

  files <- .find_ebd_sed(data_dir, ebd_file, sed_file)
  if (verbose)
    message("    清单 SED: ", basename(files$sed),
            "\n    观测 EBD: ", basename(files$ebd))

  if (verbose) message(">>> [2/8] 读取 EBD/SED 并做 zero-fill ...")
  ## 任何"读不到"的情况先返回 NULL, 交给下面的逐级自动放宽逻辑
  read_zf <- function(si, bb, yr) {
    tryCatch(.prepare_checklists(files$sed, files$ebd, species, si, yr, bb,
                                 verbose = verbose),
             error = function(e) {
               if (verbose) message("    (", conditionMessage(e), ")")
               NULL
             })
  }
  nz <- function(x) if (is.null(x)) 0 else nrow(x)
  nd <- function(x) if (is.null(x)) 0 else sum(x$species_observed)

  zf <- read_zf(sinfo, read_bb, year_range)

  ## 回退 1: 指定区域读不到 -> 放宽到全球 (说明区域与数据不匹配)
  if (nz(zf) < 30 && !auto_region) {
    zf2 <- read_zf(sinfo, global_bb, year_range)
    if (nz(zf2) > nz(zf)) {
      zf <- zf2; auto_region <- TRUE
      if (verbose) message("    指定研究区与数据位置不匹配, 改用数据实际覆盖范围。")
    }
  }
  ## 回退 2: 季节/年份样本不足 -> 放宽到全年全部年份
  if ((nz(zf) < 30 || nd(zf) < 15) && season_used != "year_round") {
    sinfo_all <- .custom_season_info("year_round", NULL, NULL, NULL, yr_place)
    zf3 <- read_zf(sinfo_all, global_bb, NULL)
    if (nd(zf3) > nd(zf) || nz(zf3) > nz(zf)) {
      zf <- zf3; sinfo <- sinfo_all; season_used <- "year_round"
      if (verbose) message("    指定时段内样本不足, 已自动改用数据全部时段 (year_round)。")
    }
  }
  if (nd(zf) < 3) {
    tops <- tryCatch(.top_species(files$ebd, 12), error = function(e) character(0))
    stop("物种 '", species, "' 在数据中检测记录不足, 无法稳定建模。\n",
         if (length(tops)) paste0("数据里记录最多的物种可能是:\n  ",
                                  paste(tops, collapse = "\n  "), "\n") else "",
         "请把 species= 改成数据中确实出现的英文名, 或换用 method=\"ebirdst\"。",
         call. = FALSE)
  }

  ## ---- 确定最终研究区 ----
  region_final <- if (auto_region) {
    .region_from_points(zf$longitude, zf$latitude)
  } else {
    pts0 <- sf::st_as_sf(zf, coords = c("longitude", "latitude"), crs = 4326,
                         remove = FALSE)
    rough <- sf::st_buffer(sf::st_transform(region_sf, 4326), 0.012)
    if (sum(lengths(sf::st_intersects(pts0, rough)) > 0) < 10) {
      if (verbose) message("    数据大多落在指定研究区之外, 改用数据实际覆盖范围。")
      .region_from_points(zf$longitude, zf$latitude)
    } else region_sf
  }
  laea_crs <- if (crs == "auto") .auto_laea(region_final) else crs
  region_laea <- sf::st_transform(region_final, laea_crs)

  ## 精确空间裁剪 (区域外扩 1 km)
  pts <- sf::st_as_sf(zf, coords = c("longitude", "latitude"), crs = 4326,
                      remove = FALSE)
  rbuf <- sf::st_transform(sf::st_buffer(region_laea, 1000), 4326)
  zf <- zf[lengths(sf::st_intersects(pts, rbuf)) > 0, , drop = FALSE]
  if (nd(zf) < 3) { ## 保险: 裁剪后太少则完全按数据范围
    region_final <- .region_from_points(pts$longitude, pts$latitude)
    laea_crs <- if (crs == "auto") .auto_laea(region_final) else crs
    region_laea <- sf::st_transform(region_final, laea_crs)
  }

  ## 预测代表日年份 -> 数据中位年份
  sinfo <- .set_rep_year(sinfo, round(stats::median(zf$year, na.rm = TRUE)))

  if (verbose)
    message(sprintf("    合格清单 %d 条, 其中检测到 %d 条 (%.1f%%)",
                    nrow(zf), sum(zf$species_observed),
                    100 * mean(zf$species_observed)))

  ## ---- train/test 划分 (小样本自动减少测试比例) ----
  zf$type <- ifelse(stats::runif(nrow(zf)) <= (1 - test_fraction), "train", "test")
  if (sum(zf$type == "train") < 60) {
    k <- max(5, floor(0.1 * nrow(zf)))
    zf$type <- "train"
    zf$type[sample(seq_len(nrow(zf)), min(k, nrow(zf)))] <- "test"
  }

  if (verbose) message(">>> [3/8] 时空网格分层子采样 ...")
  zf_ss <- .adaptive_sample(zf, laea_crs, ss_cellsize, ss_weeks, verbose)

  if (verbose) message(">>> [4/8] 在线下载环境栅格 ...")
  env <- .download_env(region_sf = region_final, laea_crs = laea_crs,
                       cache_dir = env_cache_dir, elev_z = elev_z,
                       lc_agg = lc_agg, sample_radius = sample_radius,
                       verbose = verbose)

  if (verbose) message(">>> [5/8] 提取清单点周边环境变量 ...")
  pts <- sf::st_as_sf(zf_ss, coords = c("longitude", "latitude"), crs = 4326,
                      remove = FALSE)
  env_pts <- .extract_env(pts, id_col = "checklist_id",
                          dem = env$elevation, lc = env$landcover,
                          laea_crs = laea_crs, radius = sample_radius,
                          compute_ed = compute_ed)
  zf_ss <- dplyr::left_join(zf_ss, env_pts, by = "checklist_id")

  if (verbose) message(">>> [6/8] 训练 hurdle 模型 ...")
  pland_cols <- grep("^pland_", names(zf_ss), value = TRUE)
  ed_cols    <- grep("^ed_",    names(zf_ss), value = TRUE)
  elev_cols  <- c("elevation_mean", "elevation_sd")
  effort_cols <- c("year", "day_of_year", "hours_of_day",
                   "effort_hours", "effort_distance_km",
                   "effort_speed_kmph", "number_observers")
  er_preds <- c(effort_cols, pland_cols, ed_cols, elev_cols)
  er_preds <- er_preds[er_preds %in% names(zf_ss)]

  model_df <- zf_ss %>%
    dplyr::filter(type == "train") %>%
    dplyr::select(species_observed, observation_count, all_of(er_preds)) %>%
    tidyr::drop_na(all_of(er_preds))
  if (verbose) message("    进入模型的训练清单: ", nrow(model_df), " 条")
  if (sum(model_df$species_observed) < 3)
    stop("子采样后检测记录过少, 无法建模。可增大 ss_cellsize 或检查数据。",
         call. = FALSE)

  fit <- .fit_hurdle(model_df, er_preds = er_preds, n_trees = n_trees,
                     seed = seed, hours_peak = hours_peak, verbose = verbose)

  if (verbose) message(">>> [7/8] 构建 ", grid_res, "m 预测网格并预测丰度 ...")
  pg <- .build_pred_grid(region_laea = region_laea, laea_crs = laea_crs,
                         grid_res = grid_res, dem = env$elevation,
                         lc = env$landcover, radius = sample_radius,
                         compute_ed = compute_ed,
                         rep_date = sinfo$rep_date, t_peak = fit$t_peak)
  r_layers <- .predict_to_grid(pg = pg, fit = fit, er_preds = er_preds,
                               template = pg$template)
  if (isTRUE(smooth)) {
    r_layers <- .smooth_layers(r_layers, window = 3)
    if (verbose) message("    已对丰度做轻度空间平滑 (3x3) 以降低小样本噪声。")
  }

  if (verbose) message(">>> [8/8] 出图并保存 ...")
  final_title <- if (season_used != season)
    paste0(species, " — ", .season_label(season_used, start_date, end_date, week),
           " (时段自动放宽)") else map_title
  if (do_plot) {
    .plot_rel_abundance(r_layers[["abundance"]], region_laea,
                        quantile_breaks = TRUE, n_quantiles = 10,
                        palette = "ebirdst", map_title = final_title,
                        verbose = FALSE)
  }
  .save_outputs(r_layers[["abundance"]], species, season_used, out_dir,
                out_prefix, week = week, start_date = start_date,
                end_date = end_date, verbose = verbose)
  .save_raster_layers(r_layers, species, season_used, out_dir, out_prefix,
                      week, start_date, end_date, verbose)

  diag <- .custom_diagnostics(zf_ss, fit, er_preds)
  if (verbose)
    message("    测试集遭遇率 RMSE = ", round(diag$er_rmse, 3),
            " | 计数相关 r = ", round(diag$count_cor, 3))

  structure(list(
    r = r_layers[["abundance"]],
    r_layers = r_layers,
    region = region_laea,
    species = species,
    season = .season_label(season_used, start_date, end_date, week),
    season_requested = season, season_used = season_used,
    method = "custom",
    models = fit,
    n_checklists = nrow(zf),
    n_model = nrow(model_df),
    diagnostics = diag,
    call = match.call()
  ), class = "rel_abundance")
}

## 按数据点经纬度范围构建研究区多边形 (外扩自适应, 度)
.region_from_points <- function(lon, lat, min_pad_deg = 0.05) {
  lon <- lon[is.finite(lon)]; lat <- lat[is.finite(lat)]
  bb <- sf::st_bbox(c(xmin = min(lon), xmax = max(lon),
                      ymin = min(lat), ymax = max(lat)), crs = 4326)
  pad <- max(min_pad_deg,
             0.08 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]))
  pad <- min(pad, 2)
  bb <- bb + c(-pad, -pad, pad, pad)
  sf::st_sf(name = "Data coverage", geometry = sf::st_as_sfc(bb), crs = 4326)
}

## 把季节代表日的年份替换成数据实际 (中位) 年份
.set_rep_year <- function(sinfo, yr) {
  if (!is.null(sinfo$rep_date)) {
    d <- sinfo$rep_date
    lubridate::year(d) <- as.integer(yr)
    sinfo$rep_date <- d
  }
  sinfo
}

## 小样本自适应子采样: 逐步放大时空格子, 保证训练样本量; 仍不足则用全部清单
.adaptive_sample <- function(zf, laea_crs, cellsize, weeks, verbose) {
  target_train <- 120
  cs <- cellsize; wk <- weeks
  zs <- .grid_sample_stratified(zf, laea_crs, cs, wk)
  for (i in seq_len(5)) {
    if (sum(zs$type == "train") >= target_train ||
        nrow(zs) >= 0.9 * nrow(zf)) break
    cs <- cs * 2; wk <- min(wk * 2, 13)
    zs <- .grid_sample_stratified(zf, laea_crs, cs, wk)
  }
  if (sum(zs$type == "train") < 40) {
    if (verbose) message("    样本较少, 跳过子采样, 使用全部 ", nrow(zf), " 条合格清单。")
    zs <- zf
  } else if (verbose) {
    message(sprintf("    子采样后 %d 条 (训练 %d / 测试 %d)",
                    nrow(zs), sum(zs$type == "train"), sum(zs$type == "test")))
  }
  zs
}

## 统计 EBD 中记录最多的物种 (用于物种名写错/无记录时给出友好建议)
.top_species <- function(ebd_path, n = 12) {
  cb <- readr::DataFrameCallback$new(function(x, pos) {
    names(x) <- .ebd_name_repair(names(x))
    need <- intersect(c("common_name", "category", "all_species_reported"),
                      names(x))
    x[, need, drop = FALSE]
  })
  d <- readr::read_tsv_chunked(
    ebd_path, callback = cb, chunk_size = 500000,
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE, show_col_types = FALSE, quote = "\"",
    na = c("", "NA"))
  if (nrow(d) == 0 || !"common_name" %in% names(d)) return(character(0))
  ok <- rep(TRUE, nrow(d))
  if ("category" %in% names(d)) ok <- ok & d$category == "species"
  if ("all_species_reported" %in% names(d))
    ok <- ok & .is_true(d$all_species_reported)
  tab <- sort(table(d$common_name[ok]), decreasing = TRUE)
  names(utils::head(tab, n))
}

## ----------------------------------------------------------------------------
## 3.1 时段 -> 月份集合 + 预测代表日
## ----------------------------------------------------------------------------
.custom_season_info <- function(season, week, start_date, end_date, rep_year) {
  switch(season,
    breeding = list(kind = "months", months = 5:7,
                    rep_date = as.Date(paste0(rep_year, "-06-15"))),
    nonbreeding = list(kind = "months", months = c(12, 1, 2),
                       rep_date = as.Date(paste0(rep_year, "-01-15"))),
    prebreeding_migration = list(kind = "months", months = 3:4,
                                 rep_date = as.Date(paste0(rep_year, "-04-01"))),
    postbreeding_migration = list(kind = "months", months = 8:11,
                                  rep_date = as.Date(paste0(rep_year, "-09-15"))),
    year_round = list(kind = "months", months = 1:12,
                      rep_date = as.Date(paste0(rep_year, "-06-15"))),
    weekly = {
      wk <- as.integer(week %||% 24)
      d0 <- as.Date(paste0(rep_year, "-01-01")) + (wk - 1) * 7
      list(kind = "weekly", doy1 = (wk - 1) * 7 + 1,
           doy2 = min(wk * 7, 365),
           rep_date = d0 + 3)
    },
    custom = {
      s <- as.Date(start_date); e <- as.Date(end_date)
      list(kind = "range", start = s, end = e, rep_date = s + (e - s) / 2)
    },
    stop("未知 season: ", season, call. = FALSE)
  )
}

## ----------------------------------------------------------------------------
## 3.2 在 data_dir 自动探测 EBD(观测) 与 SED(清单) 文件
## ----------------------------------------------------------------------------
.find_ebd_sed <- function(data_dir, ebd_file = NULL, sed_file = NULL) {
  if (!dir.exists(data_dir))
    stop("数据目录不存在: ", normalizePath(data_dir, mustWork = FALSE),
         "\\n请把 eBird 下载的两个 .txt 放到该目录 (可放在其任意子目录)。",
         call. = FALSE)
  ## recursive=TRUE: 兼容 data/ebd-datafile-SAMPLE/ 这类子目录
  all_txt <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$",
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)
  if (length(all_txt) == 0)
    stop("在 ", normalizePath(data_dir), " 及其子目录里没有找到 .txt 数据文件。",
         call. = FALSE)

  is_sed <- grepl("sampling", basename(all_txt), ignore.case = TRUE)
  ## EBD: 以 ebd 开头、不带 sampling; 顺带排除 BCR/IBA/USFWS/Protocol 等代码表
  is_ebd <- grepl("^ebd", basename(all_txt), ignore.case = TRUE) & !is_sed &
            !grepl("codes?\\.txt$|protocols?\\.txt$|metadata|citation|terms",
                   basename(all_txt), ignore.case = TRUE)

  sed <- if (!is.null(sed_file)) .resolve_user_file(data_dir, sed_file) else {
    if (sum(is_sed) >= 1) all_txt[is_sed][1] else NA_character_
  }
  ebd <- if (!is.null(ebd_file)) .resolve_user_file(data_dir, ebd_file) else {
    if (sum(is_ebd) >= 1) all_txt[is_ebd][1] else NA_character_
  }

  if (is.na(sed) || is.na(ebd)) {
    stop(
      "\\u65e0\\u6cd5\\u540c\\u65f6\\u627e\\u5230 \\u89c2\\u6d4b\\u6587\\u4ef6(EBD) \\u548c \\u6e05\\u5355\\u6587\\u4ef6(SED)\\u3002\\n",
      "  \\u76ee\\u5f55: ", normalizePath(data_dir), "\\n",
      "  \\u627e\\u5230\\u7684 .txt \\u6587\\u4ef6:\\n    ",
      paste(basename(all_txt), collapse = "\\n    "),
      "\\n\\n\\u8bf4\\u660e: \\u4ece\\u96f6\\u5efa\\u6a21\\u5fc5\\u987b\\u540c\\u65f6\\u6709\\u4e24\\u4e2a\\u6587\\u4ef6\\u2014\\u2014\\n",
      "  * \\u6e05\\u5355\\u6587\\u4ef6 SED: \\u6587\\u4ef6\\u540d\\u5e26 'sampling';\\n",
      "  * \\u89c2\\u6d4b\\u6587\\u4ef6 EBD: \\u6587\\u4ef6\\u540d\\u4ee5 'ebd' \\u5f00\\u5934\\u4f46\\u4e0d\\u5e26 'sampling'\\u3002\\n",
      "\\u8bf7\\u5728 eBird Custom Download \\u52fe\\u9009 'Include sampling event data'\\u3002",
      call. = FALSE)
  }
  list(ebd = ebd, sed = sed)
}

## 用户显式给的文件名: 既支持相对 data_dir, 也支持绝对路径
.resolve_user_file <- function(data_dir, f) {
  cand <- if (file.exists(f)) f else file.path(data_dir, f)
  if (!file.exists(cand)) {
    hit <- list.files(data_dir, pattern = paste0("(^|/)\\Q", basename(f), "\\E$"),
                      full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    cand <- if (length(hit)) hit[1] else cand
  }
  cand
}

## ----------------------------------------------------------------------------
## 3.3 分块读取 EBD/SED (readr, 内存安全), 边读边过滤
## ----------------------------------------------------------------------------
## 把 EBD 原始列名(大写带空格)修成小写下划线
.ebd_name_repair <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

## 判断 "完整清单" 字段 (EBD 里可能是 1/0 或 TRUE/FALSE)
.is_true <- function(v) {
  toupper(as.character(v)) %in% c("1", "TRUE", "T")
}

## 通用分块读取器
.read_ebd_chunked <- function(path, is_sampling, species, sinfo,
                              year_range, bbox_ll = NULL, chunk_size = 250000) {
  ## 只保留需要的列 (不同版本 EBD/SED 列名略有差异, 用 intersect 容错)
  keep_cols <- c("common_name", "scientific_name", "category", "observation_count",
                 "state_code", "country_code", "locality_id",
                 "latitude", "longitude", "observation_date",
                 "time_observations_started", "observer_id",
                 "sampling_event_identifier",
                 "protocol_type", "protocol_code", "protocol_name",
                 "duration_minutes", "effort_distance_km",
                 "number_observers", "all_species_reported",
                 "group_identifier")

  cb <- readr::DataFrameCallback$new(function(x, pos) {
    names(x) <- .ebd_name_repair(names(x))
    x <- x[, intersect(keep_cols, names(x)), drop = FALSE]
    if (nrow(x) == 0) return(x[0, , drop = FALSE])

    ## 只要完整清单 (all_species_reported)
    x <- x[.is_true(x$all_species_reported), , drop = FALSE]
    if (nrow(x) == 0) return(x[0, , drop = FALSE])

    x$latitude  <- as.numeric(x$latitude)
    x$longitude <- as.numeric(x$longitude)
    x <- x[!is.na(x$latitude) & !is.na(x$longitude), , drop = FALSE]

    ## 可选空间裁剪: 仅当显式传入 bbox_ll 时才做 (region="auto" 时边界未知, 不裁)
    if (!is.null(bbox_ll)) {
      x <- x[x$longitude >= bbox_ll["xmin"] & x$longitude <= bbox_ll["xmax"] &
             x$latitude  >= bbox_ll["ymin"] & x$latitude  <= bbox_ll["ymax"], ,
             drop = FALSE]
      if (nrow(x) == 0) return(x[0, , drop = FALSE])
    }

    x$observation_date <- as.Date(x$observation_date)
    x <- x[!is.na(x$observation_date), , drop = FALSE]
    yr <- as.integer(format(x$observation_date, "%Y"))
    mo <- as.integer(format(x$observation_date, "%m"))
    doy <- as.integer(format(x$observation_date, "%j"))
    ## year_range=NULL 表示不限制年份 (自动适配数据实际年份)
    yr_keep <- if (is.null(year_range)) rep(TRUE, nrow(x)) else yr %in% year_range
    keep_date <- switch(sinfo$kind,
      months = yr_keep & mo %in% sinfo$months,
      weekly = yr_keep & doy >= sinfo$doy1 & doy <= sinfo$doy2,
      range  = x$observation_date >= sinfo$start &
               x$observation_date <= sinfo$end,
      rep(TRUE, nrow(x)))
    x <- x[keep_date, , drop = FALSE]
    if (nrow(x) == 0) return(x[0, , drop = FALSE])

    ## EBD: 只保留目标物种 (且 category 为 species; 排除 spuh/杂交/家型)
    if (!is_sampling && !is.null(species)) {
      ## 同时接受通用名(Wood Thrush)和学名(Hylocichla mustelina)
      sp_lc  <- tolower(trimws(species))
      hit_cn <- tolower(trimws(x$common_name)) == sp_lc
      hit_sn <- if ("scientific_name" %in% names(x))
        tolower(trimws(x$scientific_name)) == sp_lc else rep(FALSE, nrow(x))
      x <- x[(hit_cn | hit_sn) &
             (x$category %in% c("species") | is.na(x$category)), ,
           drop = FALSE]
    }
    x
  })

  readr::read_tsv_chunked(
    path, callback = cb,
    col_types = readr::cols(.default = readr::col_character()),
    chunk_size = chunk_size, progress = FALSE, show_col_types = FALSE,
    quote = "\"", na = c("", "NA"))
}

## 开始时刻 "HH:MM:SS" -> 小数小时
.time_to_decimal <- function(x) {
  tt <- suppressWarnings(hms::as_hms(x))
  lubridate::hour(tt) + lubridate::minute(tt) / 60 +
    lubridate::second(tt) / 3600
}

## zero-fill + 努力量派生 + effort 过滤 (空间精确裁剪在 .run_custom 中完成)
.prepare_checklists <- function(sed_path, ebd_path, species, sinfo,
                                year_range, bbox_ll, verbose) {
  ## 读清单(SED)
  chk <- .read_ebd_chunked(sed_path, is_sampling = TRUE, species = NULL,
                           sinfo = sinfo, year_range = year_range,
                           bbox_ll = bbox_ll)
  if (nrow(chk) == 0) stop("SED 清单文件在该时段/地区没有任何记录。", call. = FALSE)
  chk <- chk[!duplicated(chk$sampling_event_identifier), , drop = FALSE]
  chk <- dplyr::rename(chk, checklist_id = sampling_event_identifier)

  ## 读观测(EBD)
  obs <- .read_ebd_chunked(ebd_path, is_sampling = FALSE, species = species,
                           sinfo = sinfo, year_range = year_range,
                           bbox_ll = bbox_ll)
  if (nrow(obs) == 0)
    stop("EBD 中找不到物种 '", species, "' 在该时段/地区的记录。", call. = FALSE)

  obs$count_num <- suppressWarnings(as.integer(obs$observation_count))
  counts <- obs %>%
    dplyr::group_by(sampling_event_identifier) %>%
    dplyr::summarise(
      seen = TRUE,
      count_val = ifelse(all(is.na(count_num)), NA_integer_,
                         as.integer(sum(count_num, na.rm = TRUE))),
      .groups = "drop") %>%
    dplyr::rename(checklist_id = sampling_event_identifier)

  zf <- dplyr::left_join(chk, counts, by = "checklist_id")
  zf$species_observed <- !is.na(zf$seen)
  zf$seen <- NULL
  zf$observation_count <- zf$count_val
  zf$count_val <- NULL

  ## ---- 协议字段容错 -------------------------------------------------------
  ## 教程只保留 Stationary(P20) / Traveling(P21); 但真实 SAMPLE 数据常只有
  ## PROTOCOL CODE, 且大量是 Exhaustive Area(P22, 区域穷举计数)。工程上把
  ## P22/P23 归入"行进类"纳入 (它们同样有 duration, 是高质量存在/缺失数据),
  ## 缺失距离再用中位数兜底, 避免小样本被整个删光。
  if (!"protocol_type" %in% names(zf)) {
    if ("protocol_code" %in% names(zf)) {
      pmap <- c(P20 = "Stationary", P21 = "Traveling",
                P22 = "Traveling", P23 = "Traveling")
      pt <- unname(pmap[trimws(as.character(zf$protocol_code))])
      pt[is.na(pt)] <- "Traveling"   # 未知协议: 有时长就按行进类纳入
      zf$protocol_type <- pt
      if (isTRUE(verbose))
        message("    (数据无 PROTOCOL TYPE, 已按 PROTOCOL CODE 映射; P22 区域穷举归入行进类)")
    } else {
      zf$protocol_type <- "Traveling"
    }
  } else {
    zf$protocol_type[is.na(zf$protocol_type)] <- "Traveling"
  }

  zf$duration_minutes   <- as.numeric(zf$duration_minutes)
  zf$effort_distance_km <- as.numeric(zf$effort_distance_km)
  zf$number_observers   <- as.integer(zf$number_observers)
  zf$latitude           <- as.numeric(zf$latitude)
  zf$longitude          <- as.numeric(zf$longitude)

  ## 定点调查距离记 0; 行进类距离缺失则用中位数兜底 (全缺省 1 km), 不再删行
  zf$effort_distance_km[zf$protocol_type == "Stationary"] <- 0
  is_trav <- zf$protocol_type == "Traveling"
  dmed <- stats::median(zf$effort_distance_km[is_trav], na.rm = TRUE)
  if (!is.finite(dmed)) dmed <- 1
  zf$effort_distance_km[is_trav & is.na(zf$effort_distance_km)] <- dmed

  zf <- zf %>% dplyr::mutate(
    effort_hours       = duration_minutes / 60,
    effort_speed_kmph  = ifelse(is.na(effort_hours) | effort_hours == 0,
                                NA_real_, effort_distance_km / effort_hours),
    hours_of_day       = .time_to_decimal(time_observations_started),
    year               = as.integer(format(observation_date, "%Y")),
    day_of_year        = as.integer(format(observation_date, "%j"))
  )

  ## 速度兜底 (理论上距离兜底后已无 NA; 双保险)
  smed <- stats::median(zf$effort_speed_kmph, na.rm = TRUE)
  if (!is.finite(smed)) smed <- 2
  zf$effort_speed_kmph[is.na(zf$effort_speed_kmph)] <- smed
  ## 开始时刻兜底
  if (all(is.na(zf$hours_of_day))) zf$hours_of_day <- 7
  med_h <- stats::median(zf$hours_of_day, na.rm = TRUE)
  if (!is.finite(med_h)) med_h <- 7
  zf$hours_of_day[is.na(zf$hours_of_day)] <- med_h
  zf$number_observers[is.na(zf$number_observers)] <- 1L

  ## 只保留硬性质量条件: 协议类型、时长 0-6h、距离/速度合理、坐标与人数;
  ## 距离已兜底, 不再因缺距离删行
  zf <- zf %>% dplyr::filter(
    protocol_type %in% c("Stationary", "Traveling"),
    !is.na(effort_hours), effort_hours > 0, effort_hours <= 6,
    is.na(effort_distance_km) | effort_distance_km <= 10,
    is.na(effort_speed_kmph)  | effort_speed_kmph  <= 100,
    is.na(number_observers)   | number_observers   <= 10,
    !is.na(latitude), !is.na(longitude)
  )
  zf
}


## ----------------------------------------------------------------------------
## 3.4 时空网格分层子采样 (复刻 ebirdst::grid_sample_stratified 的思想)
##     在 每个 3km×3km×1周×年 的格子里, 检测/非检测各随机抽 1 条
## ----------------------------------------------------------------------------
.grid_sample_stratified <- function(df, laea_crs, cellsize = 3000, weeks = 1) {
  pts <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326,
                  remove = FALSE) %>% st_transform(laea_crs)
  cc <- st_coordinates(pts)
  df$.gx <- floor(cc[, 1] / cellsize)
  df$.gy <- floor(cc[, 2] / cellsize)
  df$.gw <- floor((df$day_of_year - 1) / (7 * weeks))
  df %>%
    group_by(type, species_observed, .gx, .gy, year, .gw) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-.gx, -.gy, -.gw)
}

## ----------------------------------------------------------------------------
## 3.5 在线下载环境栅格 (免密钥)
##     高程: elevatr (AWS Terrain Tiles, 约 1km, 对齐教程 GMTED)
##     土地覆盖: ESA WorldCover 2021 10m (Planetary Computer, /vsicurl 窗口读取)
##     两个都投影到研究区 LAEA, 土地覆盖最近邻聚合到 lc_agg 米
## ----------------------------------------------------------------------------
.WC_CLASSES <- data.frame(
  code = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100),
  name = c("tree", "shrub", "grass", "cropland", "built", "bare",
           "snow", "water", "wetland", "mangrove", "moss"),
  stringsAsFactors = FALSE)

## 给 WorldCover asset 的 href 申请匿名 SAS 签名
.wc_sign <- function(href) {
  u <- paste0("https://planetarycomputer.microsoft.com/api/sas/v1/sign?href=",
              utils::URLencode(href, reserved = TRUE))
  httr::content(httr::GET(u, httr::timeout(60)), as = "parsed")$href
}

.download_env <- function(region_sf, laea_crs, cache_dir, elev_z, lc_agg,
                          sample_radius, verbose) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  bbll <- sf::st_bbox(sf::st_transform(region_sf, 4326))
  key <- sprintf("z%s_lc%s_%s", elev_z, lc_agg,
                 paste(round(bbll, 1), collapse = "_"))
  key <- gsub("[^0-9A-Za-z_.-]+", "", key)
  f_elev <- file.path(cache_dir, paste0("elevation_", key, ".tif"))
  f_lc   <- file.path(cache_dir, paste0("landcover_", key, ".tif"))

  ## ---------- 高程 (失败不致命, 后面用 0 兜底) ----------
  dem <- NULL
  if (file.exists(f_elev)) {
    dem <- terra::rast(f_elev)
    if (verbose) message("    高程(缓存): ", basename(f_elev))
  } else {
    dem <- tryCatch({
      if (verbose) message("    下载高程 elevatr (z=", elev_z, ") ...")
      region_ll <- sf::st_transform(region_sf, 4326)
      dem_r <- elevatr::get_elev_raster(locations = region_ll, z = elev_z,
                                        clip = "locations", verbose = FALSE)
      d <- terra::rast(dem_r)
      if (is.na(terra::crs(d)) || terra::crs(d) == "") terra::crs(d) <- "EPSG:4326"
      d <- terra::project(d, laea_crs, method = "bilinear", res = 1000)
      terra::writeRaster(d, f_elev, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
      d
    }, error = function(e) {
      warning("高程数据下载失败, 本次不使用高程特征: ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
  }

  ## ---------- 土地覆盖 (失败不致命, 后面用 0 兜底) ----------
  lc <- NULL
  if (file.exists(f_lc)) {
    lc <- terra::rast(f_lc)
    if (verbose) message("    土地覆盖(缓存): ", basename(f_lc))
  } else {
    lc <- tryCatch({
      if (verbose) message("    下载 ESA WorldCover 土地覆盖 (按 tile 窗口读取) ...")
      bb_poly <- sf::st_as_sfc(sf::st_bbox(bbll + c(-.05, -.05, .05, .05)),
                               crs = 4326)
      surl <- sprintf(
        "https://planetarycomputer.microsoft.com/api/stac/v1/search?collections=esa-worldcover&bbox=%s&limit=100",
        paste(sf::st_bbox(bb_poly)[c("xmin", "ymin", "xmax", "ymax")],
              collapse = ","))
      js <- jsonlite::fromJSON(
        httr::content(httr::GET(surl, httr::timeout(60)), as = "text",
                      encoding = "UTF-8"),
        simplifyVector = FALSE)
      if (length(js$features) == 0)
        stop("Planetary Computer 没有返回该区域的 WorldCover 影像。", call. = FALSE)
      if (verbose) message("      相交 tile 数: ", length(js$features))

      templ <- terra::rast(terra::vect(sf::st_transform(bb_poly, laea_crs)),
                           res = lc_agg, crs = laea_crs)
      out <- NULL
      for (i in seq_along(js$features)) {
        href <- js$features[[i]]$assets$map$href
        signed <- .wc_sign(href)
        r <- terra::rast(paste0("/vsicurl/", signed))
        win <- terra::crop(r, terra::vect(bb_poly), snap = "out")
        lc_i <- terra::project(win, templ, method = "near")
        out <- if (is.null(out)) lc_i else terra::cover(out, lc_i)
        if (verbose) message("      tile ", i, "/", length(js$features), " 完成")
      }
      out <- terra::classify(out,
        cbind(-Inf, .WC_CLASSES$code[1] - 1, NA), include.lowest = TRUE)
      terra::writeRaster(out, f_lc, overwrite = TRUE,
                         datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
      out
    }, error = function(e) {
      warning("土地覆盖数据下载失败, 本次不使用土地覆盖特征: ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
  }

  list(elevation = dem, landcover = lc)
}

## ----------------------------------------------------------------------------
## 3.6 点周边环境变量提取 (exactextractr, C++ 加速)
##     高程: mean / sd ; 土地覆盖: 各类面积百分比 pland_*
## ----------------------------------------------------------------------------
.extract_env <- function(pts_sf, id_col, dem, lc, laea_crs, radius = 1500,
                         compute_ed = FALSE) {
  pts_laea <- sf::st_transform(pts_sf, laea_crs)
  rownames(pts_laea) <- NULL
  buf <- sf::st_buffer(pts_laea, dist = radius)
  n <- nrow(pts_laea)

  ## ---- 高程 mean / sd (dem 缺失则全 0) ----
  if (is.null(dem)) {
    de <- data.frame(elevation_mean = rep(0, n), elevation_sd = rep(0, n))
  } else {
    de <- exactextractr::exact_extract(dem, buf, fun = c("mean", "stdev"),
                                       progress = FALSE)
    names(de) <- c("elevation_mean", "elevation_sd")
    de$elevation_sd[is.nan(de$elevation_sd)] <- 0
    de$elevation_mean[is.na(de$elevation_mean)] <- 0
    de$elevation_sd[is.na(de$elevation_sd)] <- 0
  }

  ## ---- 土地覆盖各类面积比例 (lc 缺失则全 0) ----
  pland <- matrix(0, nrow = n, nrow(.WC_CLASSES))
  colnames(pland) <- paste0("pland_", .WC_CLASSES$name)
  if (!is.null(lc)) {
    frac <- exactextractr::exact_extract(
      lc, buf, progress = FALSE,
      fun = function(values, coverage_fractions) {
        ok <- !is.na(values)
        if (!any(ok)) return(numeric(0))
        tab <- tapply(coverage_fractions[ok], values[ok], sum, na.rm = TRUE)
        tab / sum(tab, na.rm = TRUE)
      })
    for (i in seq_along(frac)) {
      if (length(frac[[i]]) > 0) {
        codes <- as.integer(names(frac[[i]]))
        idx <- match(codes, .WC_CLASSES$code)
        pland[i, idx[!is.na(idx)]] <- as.numeric(frac[[i]])[!is.na(idx)]
      }
    }
  }

  out <- data.frame(id = pts_laea[[id_col]],
                    elevation_mean = de$elevation_mean,
                    elevation_sd = de$elevation_sd,
                    pland, stringsAsFactors = FALSE,
                    check.names = FALSE)
  names(out)[1] <- id_col
  out
}

## ----------------------------------------------------------------------------
## 3.7 hurdle 模型训练
## ----------------------------------------------------------------------------
## MCC-F1 最佳阈值; 兼容 mccf1 旧版与 1.x 新版返回结构
.best_mccf1_threshold <- function(response, predictor) {
  m <- mccf1::mccf1(response = response, predictor = predictor)
  ## 旧版本 (<1.0): list(mccf1_summary = data.frame(best_threshold=...))
  if (!is.null(m$mccf1_summary) &&
      "best_threshold" %in% names(m$mccf1_summary)) {
    return(m$mccf1_summary$best_threshold[1])
  }
  ## 新版本 (1.x): 三个等长向量 normalized_mcc / f1 / thresholds
  score <- m$normalized_mcc
  thr <- m$thresholds
  ok <- is.finite(score) & is.finite(thr)
  if (!any(ok)) return(0.5)
  thr[ok][which.max(score[ok])]
}

.fit_hurdle <- function(model_df, er_preds, n_trees, seed,
                        hours_peak, verbose) {
  er_df <- model_df %>% dplyr::select(species_observed, all_of(er_preds))
  n1 <- sum(er_df$species_observed)
  n0 <- sum(!er_df$species_observed)
  detection_freq <- mean(er_df$species_observed)

  er_model <- NULL; threshold <- 0; cal_model <- NULL
  if (n0 == 0) {
    ## 全部检测到 (极罕见): 遭遇率恒为 1
    if (verbose) message("    训练集中几乎全部为检测记录, 遭遇率按常数处理。")
  } else {
    er_model <- ranger::ranger(
      as.factor(species_observed) ~ ., data = er_df,
      probability = TRUE, replace = TRUE,
      sample.fraction = c(detection_freq, detection_freq),
      importance = "impurity", num.trees = n_trees, seed = seed,
      verbose = FALSE)
    oob_er <- er_model$predictions[, 2]
    ## 小样本 + 高检测率时, 少数样本可能从未成为 OOB -> 用模型自预测填补
    if (anyNA(oob_er)) {
      ins <- tryCatch(predict(er_model, data = er_df)$predictions[, 2],
                      error = function(e) rep(detection_freq, nrow(er_df)))
      oob_er[is.na(oob_er)] <- ins[is.na(oob_er)]
      oob_er[is.na(oob_er)] <- detection_freq
    }
    threshold <- .best_mccf1_threshold(as.integer(er_df$species_observed), oob_er)
    ## SCAM 校准; 失败则退回原始概率
    cal_model <- tryCatch(
      scam::scam(obs ~ s(pred, k = 6, bs = "mpi"),
                 data = data.frame(obs = as.integer(er_df$species_observed),
                                   pred = oob_er), gamma = 2),
      error = function(e) {
        warning("SCAM 校准失败, 使用未校准概率: ", conditionMessage(e),
                call. = FALSE)
        NULL
      })
  }

  ## ---------- 第二关: 计数随机森林 ----------
  count_model <- NULL; count_const <- NA_real_
  if (!is.null(er_model)) {
    count_er <- predict(er_model, data = er_df)$predictions[, 2]
    count_er[is.na(count_er)] <- detection_freq
    cnt_idx <- which(!is.na(model_df$observation_count) &
                     (model_df$observation_count > 0 | count_er > threshold))
  } else {
    count_er <- rep(1, nrow(er_df))
    cnt_idx <- which(!is.na(model_df$observation_count) &
                     model_df$observation_count > 0)
  }
  count_df <- model_df[cnt_idx, , drop = FALSE]
  count_preds <- c(er_preds, "predicted_er")
  if (!is.null(er_model)) count_df$predicted_er <- count_er[cnt_idx]
  count_train <- count_df %>%
    dplyr::select(observation_count, all_of(intersect(count_preds, names(count_df)))) %>%
    tidyr::drop_na()

  if (nrow(count_train) >= 10 &&
      length(unique(count_train$observation_count)) >= 2) {
    count_model <- ranger::ranger(
      observation_count ~ ., data = count_train,
      replace = TRUE, importance = "impurity",
      num.trees = n_trees, seed = seed, verbose = FALSE)
  } else {
    pos <- model_df$observation_count
    pos <- pos[!is.na(pos) & pos > 0]
    count_const <- if (length(pos)) mean(pos) else 1
    if (verbose) message("    计数样本不足, 计数层使用常数均值 ",
                         round(count_const, 2))
  }

  t_peak <- if (!is.null(hours_peak)) {
    as.numeric(hours_peak)
  } else if (!is.null(er_model)) {
    tryCatch(.estimate_peak_hours(er_model, cal_model, er_df),
             error = function(e) 7)
  } else 7
  if (verbose) message("    最佳观测开始时刻 hours_of_day = ", round(t_peak, 2),
                       " ; 遭遇率阈值 = ", round(threshold, 3))

  list(er_model = er_model, cal_model = cal_model,
       count_model = count_model, count_const = count_const,
       threshold = threshold, count_preds = count_preds,
       er_preds = er_preds, t_peak = t_peak,
       detection_freq = detection_freq)
}

## ----------------------------------------------------------------------------
## 3.8 预测网格: LAEA 3km 模板 + 中心点环境 + 标准努力量
## ----------------------------------------------------------------------------
.build_pred_grid <- function(region_laea, laea_crs, grid_res, dem, lc,
                             radius, compute_ed, rep_date, t_peak) {
  ## 栅格模板, 区内填 1
  template <- terra::rast(region_laea, res = rep(grid_res, 2))
  template <- terra::rasterize(terra::vect(region_laea), template, field = 1)
  names(template) <- "study_region"

  cells <- as.data.frame(template, cells = TRUE, xy = TRUE, na.rm = TRUE)
  names(cells) <- c("cell", "x", "y", "study_region")

  pts <- st_as_sf(cells, coords = c("x", "y"), crs = laea_crs, remove = FALSE)
  env_pg <- .extract_env(pts, id_col = "cell", dem = dem, lc = lc,
                         laea_crs = laea_crs, radius = radius,
                         compute_ed = compute_ed)
  g <- left_join(cells, env_pg, by = "cell")

  ## 标准努力量 (复刻教程): 1 小时、2 km、2 km/h、1 名观察者、最佳时刻
  g <- g %>% mutate(
    observation_date = rep_date,
    year = as.integer(format(rep_date, "%Y")),
    day_of_year = as.integer(format(rep_date, "%j")),
    hours_of_day = t_peak,
    effort_distance_km = 2,
    effort_hours = 1,
    effort_speed_kmph = 2,
    number_observers = 1)
  list(grid = g, template = template)
}

## ----------------------------------------------------------------------------
## 3.9 在网格上预测, 栅格化四层结果
## ----------------------------------------------------------------------------
.predict_to_grid <- function(pg, fit, er_preds, template) {
  g <- pg$grid
  need <- er_preds
  for (cc in setdiff(need, names(g))) g[[cc]] <- 0

  ## 遭遇率原始概率 (er_model 为空时恒为 1)
  if (is.null(fit$er_model)) {
    raw <- rep(1, nrow(g))
  } else {
    raw <- predict(fit$er_model, data = g[, need, drop = FALSE])$predictions[, 2]
    raw[is.na(raw)] <- if (!is.null(fit$detection_freq)) fit$detection_freq else 0
  }
  in_range <- as.integer(raw > fit$threshold)

  ## 校准 (cal_model 为空时用原始概率), 截断 0-1
  if (is.null(fit$cal_model)) {
    cal <- raw
  } else {
    cal <- as.numeric(predict(fit$cal_model,
                              newdata = data.frame(pred = raw),
                              type = "response"))
  }
  cal[is.na(cal)] <- raw[is.na(cal)]
  cal[cal < 0] <- 0; cal[cal > 1] <- 1

  ## 计数
  g$predicted_er <- raw
  if (is.null(fit$count_model)) {
    cnt <- rep(fit$count_const, nrow(g))
  } else {
    cneed <- fit$count_preds
    for (cc in setdiff(cneed, names(g))) g[[cc]] <- 0
    cnt <- predict(fit$count_model,
                   data = g[, cneed, drop = FALSE])$predictions
  }
  fallback_cnt <- if (is.finite(fit$count_const)) fit$count_const else 1
  cnt[is.na(cnt) | cnt < 0] <- fallback_cnt

  abundance <- in_range * cal * cnt

  pred <- data.frame(cell = g$cell, x = g$x, y = g$y,
                     in_range = in_range, encounter_rate = cal,
                     count = cnt, abundance = abundance)
  pts <- sf::st_as_sf(pred, coords = c("x", "y"), crs = terra::crs(template),
                      remove = FALSE)
  r <- terra::rasterize(terra::vect(pts), template,
                        field = c("in_range", "encounter_rate",
                                  "count", "abundance"))
  r
}

## ============================================================================
## 4. 诊断 / 研究区解析 / 投影 / 绘图 / 保存 / S3
## ============================================================================

## 简单测试集诊断 (对 er_model/cal_model/count_model 为 NULL 的兜底情形安全)
.custom_diagnostics <- function(zf_ss, fit, er_preds) {
  er_rmse <- NA_real_; count_cor <- NA_real_
  test <- zf_ss %>% dplyr::filter(type == "test")
  test <- tidyr::drop_na(test, tidyselect::all_of(intersect(er_preds, names(test))))
  if (nrow(test) >= 5) {
    if (is.null(fit$er_model)) {
      cal <- rep(1, nrow(test))
    } else {
      ep <- intersect(er_preds, names(test))
      raw <- tryCatch(
        predict(fit$er_model, data = test[, ep, drop = FALSE])$predictions[, 2],
        error = function(e) rep(NA_real_, nrow(test)))
      if (is.null(fit$cal_model)) cal <- raw else
        cal <- tryCatch(
          as.numeric(predict(fit$cal_model, newdata = data.frame(pred = raw),
                             type = "response")),
          error = function(e) raw)
    }
    cal[is.na(cal)] <- mean(cal, na.rm = TRUE)
    er_rmse <- sqrt(mean((as.integer(test$species_observed) - cal) ^ 2,
                         na.rm = TRUE))
    det <- test[!is.na(test$observation_count) & test$observation_count > 0, ]
    if (!is.null(fit$count_model) && nrow(det) >= 10) {
      cp <- intersect(fit$count_preds, names(det))
      pc <- tryCatch(
        predict(fit$count_model, data = det[, cp, drop = FALSE])$predictions,
        error = function(e) rep(NA_real_, nrow(det)))
      ok <- is.finite(pc) & is.finite(det$observation_count)
      if (sum(ok) >= 10 && stats::sd(pc[ok]) > 0)
        count_cor <- suppressWarnings(
          stats::cor(det$observation_count[ok], pc[ok], method = "spearman"))
    }
  }
  list(er_rmse = er_rmse, count_cor = count_cor)
}

## ----------------------------------------------------------------------------
## 研究区解析 (rnaturalearth, 自动联网下载边界)
## ----------------------------------------------------------------------------
.resolve_region_ne <- function(region) {
  .check_pkgs(c("rnaturalearth", "rnaturalearthdata", "sf"))

  if (inherits(region, c("sf", "sfc"))) {
    if (!inherits(region, "sf")) region <- sf::st_sf(geometry = region)
    out <- sf::st_transform(region, 4326)
    return(out)
  }

  type <- tolower(region$type %||% "state")

  if (type %in% c("none", "auto", "data")) {
    ## custom 法: 暂用全球框占位并打标记, 读完数据后按数据实际范围定界
    g <- sf::st_sfc(sf::st_polygon(list(matrix(
      c(-180, -85, 180, -85, 180, 85, -180, 85, -180, -85),
      ncol = 2, byrow = TRUE))), crs = 4326)
    out <- sf::st_sf(name = "Data-defined", geometry = g)
    attr(out, "auto_region") <- TRUE
    return(out)
  }

  if (type == "conus") return(.conus_sf())

  if (type == "state") {
    country_iso <- region$country_iso %||% "US"
    st <- .fetch_states_sf(country_iso)
    hit <- grepl(region$name, st$name, ignore.case = TRUE)
    if (!any(hit))
      stop("找不到州/省: ", region$name, " (", country_iso, ")。",
           "可用示例: 'Georgia', 'Michigan', 'California'。", call. = FALSE)
    out <- st[hit, ]
    if (nrow(out) > 1) out <- sf::st_union(out) |> sf::st_sf()
    return(sf::st_transform(out, 4326))
  }

  if (type == "country") {
    ct <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    hit <- grepl(region$name, ct$name, ignore.case = TRUE) |
           grepl(region$name, ct$name_long, ignore.case = TRUE) |
           tolower(ct$iso_a2) == tolower(region$name) |
           tolower(ct$iso_a3) == tolower(region$name)
    if (!any(hit))
      stop("找不到国家: ", region$name,
           "。提示: 美国请写 'United States of America'。", call. = FALSE)
    out <- ct[hit, ]
    if (nrow(out) > 1) out <- sf::st_union(out) |> sf::st_sf()
    return(sf::st_transform(out, 4326))
  }

  stop("region$type 只支持 'state'/'country'/'conus'/'auto'/'none' 或 sf 对象。",
       call. = FALSE)
}

## 高精度州界需要 rnaturalearthhires (已从 CRAN archive, 走 r-universe)
.ensure_rnaturalearthhires <- function() {
  if (requireNamespace("rnaturalearthhires", quietly = TRUE))
    return(invisible(TRUE))
  message(">>> 首次使用州界, 尝试安装 rnaturalearthhires ...")
  tryCatch(utils::install.packages(
    "rnaturalearthhires",
    repos = c("https://ropensci.r-universe.dev/",
              "https://cloud.r-project.org"), quiet = TRUE),
    error = function(e) NULL)
  invisible(requireNamespace("rnaturalearthhires", quietly = TRUE))
}

## 取州/省边界; hires 不可用时回退到内置低分辨率 states
.fetch_states_sf <- function(country_iso = "US") {
  if (.ensure_rnaturalearthhires()) {
    return(rnaturalearth::ne_states(iso_a2 = country_iso,
                                    returnclass = "sf"))
  }
  message("    rnaturalearthhires 不可用, 改用低分辨率州界 ...")
  rnaturalearth::ne_download(scale = 110, type = "states",
                             category = "cultural", returnclass = "sf")
}

## 美国本土 48 州 (排除 Alaska / Hawaii / Puerto Rico)
.conus_sf <- function() {
  st <- .fetch_states_sf("US")
  keep <- !grepl("alaska|hawaii|puerto rico", st$name, ignore.case = TRUE)
  out <- sf::st_union(st[keep, ]) |> sf::st_sf()
  sf::st_transform(out, 4326)
}

## 按研究区质心自动选 Lambert 方位等积投影 (LAEA)
.auto_laea <- function(region_sf) {
  geom <- region_sf |>
    sf::st_transform(4326) |>
    sf::st_geometry() |>
    sf::st_make_valid()
  cent <- geom |> sf::st_centroid() |> sf::st_coordinates()
  lon <- mean(cent[, 1], na.rm = TRUE)
  lat <- mean(cent[, 2], na.rm = TRUE)
  paste0("+proj=laea +lat_0=", round(lat, 3),
         " +lon_0=", round(lon, 3),
         " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
}

## ----------------------------------------------------------------------------
## 丰度配色: 优先 ebirdst 官方调色板, 没装包则用近似的内置色带
## ----------------------------------------------------------------------------
.abundance_pal <- function(n) {
  if (requireNamespace("ebirdst", quietly = TRUE))
    return(ebirdst::ebirdst_palettes(n))
  grDevices::colorRampPalette(
    c("#3d0963", "#4575b4", "#1a9850", "#fee08b", "#f46d43", "#a50026"))(n)
}

##' 出图: 复刻 ebirdst applications vignette 的标准样式
.plot_rel_abundance <- function(r, region_sf,
                                quantile_breaks, n_quantiles,
                                palette, map_title, verbose) {
  if (verbose) message(">>> [4/5] 计算分箱 ...")
  vals <- terra::values(r, na.rm = TRUE, mat = FALSE)
  vals_pos <- vals[vals > 0 & !is.na(vals)]
  if (length(vals_pos) < 10) {
    warning("研究区内非零栅格单元太少, 地图可能没有意义。")
    brks <- seq(0, max(vals, na.rm = TRUE), length.out = n_quantiles + 1)
  } else if (quantile_breaks) {
    qs <- stats::quantile(vals_pos, seq(0, 1, length.out = n_quantiles + 1),
                          na.rm = TRUE)
    brks <- unique(c(0, as.numeric(qs)))
  } else {
    brks <- seq(0, max(vals, na.rm = TRUE), length.out = n_quantiles + 1)
  }
  if (identical(palette, "ebirdst")) {
    pal <- c("#e6e6e6", .abundance_pal(max(1, length(brks) - 2)))
  } else pal <- palette

  if (verbose) message(">>> [5/5] 绘图 ...")
  crs_map <- terra::crs(r)
  countries <- rnaturalearth::ne_countries(returnclass = "sf") |>
    sf::st_geometry() |> sf::st_transform(crs_map)
  countries_v <- terra::vect(countries)
  states <- tryCatch(.fetch_states_sf("US") |> sf::st_transform(crs_map),
                     error = function(e) NULL)
  states_v <- if (!is.null(states)) terra::vect(states) else NULL
  region_v <- if (!is.null(region_sf))
    terra::vect(sf::st_transform(region_sf, crs_map)) else NULL

  graphics::par(mar = c(4, 0.5, 1, 0.5))
  if (!is.null(region_v)) {
    terra::plot(region_v, col = NA, border = NA, axes = FALSE)
    terra::plot(countries_v, col = "#cfcfcf", border = "#888888",
                lwd = 0.5, add = TRUE)
  } else {
    terra::plot(countries_v, col = "#cfcfcf", border = "#888888",
                lwd = 0.5, axes = FALSE)
  }
  terra::plot(r, col = pal, breaks = brks, maxpixels = terra::ncell(r),
              legend = FALSE, axes = FALSE, add = TRUE)
  if (!is.null(states_v))
    terra::lines(states_v, col = "#ffffff", lwd = 0.75)
  terra::lines(countries_v, col = "#ffffff", lwd = 1.5)
  if (!is.null(region_v))
    terra::lines(region_v, col = "#000000", lwd = 1.2)

  labels_q <- stats::quantile(brks, c(0, 0.5, 1))
  label_breaks <- seq(0, 1, length.out = length(brks))
  fields::image.plot(
    zlim = c(0, 1), breaks = label_breaks, col = pal,
    smallplot = c(0.90, 0.93, 0.15, 0.85), legend.only = TRUE,
    axis.args = list(at = c(0, 0.5, 1), labels = round(labels_q, 2),
                     col.axis = "black", fg = NA, cex.axis = 0.9,
                     lwd.ticks = 0, line = -0.5))
  graphics::mtext(map_title, side = 3, line = -1, cex = 1.1)
  invisible(NULL)
}

## ----------------------------------------------------------------------------
## 季节英文标签 (用于标题与文件名)
## ----------------------------------------------------------------------------
.season_label <- function(season, start_date, end_date, week = NULL) {
  switch(season,
    breeding               = "Breeding Season",
    nonbreeding            = "Non-breeding Season",
    prebreeding_migration  = "Pre-breeding Migration",
    postbreeding_migration = "Post-breeding Migration",
    year_round             = "Full Year",
    weekly                 = paste0("Week ", week),
    custom                 = paste0(start_date, "_to_", end_date),
    season)
}

## 生成文件系统安全的输出前缀
.safe_prefix <- function(species, season, week = NULL,
                         start_date = NULL, end_date = NULL) {
  lab <- .season_label(season, start_date, end_date, week)
  p <- paste(gsub("[^A-Za-z0-9]+", "-", species),
             gsub("[^A-Za-z0-9]+", "-", lab), sep = "-")
  gsub("^-|-$", "", gsub("-+", "-", p))
}

## 保存单层丰度 GeoTIFF
.save_outputs <- function(r, species, season, out_dir, out_prefix = NULL,
                          week = NULL, start_date = NULL, end_date = NULL,
                          verbose = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(out_prefix))
    out_prefix <- .safe_prefix(species, season, week, start_date, end_date)
  f <- file.path(out_dir, paste0(out_prefix, "_abundance.tif"))
  terra::writeRaster(r, f, overwrite = TRUE)
  if (verbose) message("    已保存丰度栅格: ", f)
  invisible(f)
}

## 保存四层 (in_range/encounter_rate/count/abundance) GeoTIFF
.save_raster_layers <- function(r_layers, species, season, out_dir,
                                out_prefix = NULL, week = NULL,
                                start_date = NULL, end_date = NULL,
                                verbose = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(out_prefix))
    out_prefix <- .safe_prefix(species, season, week, start_date, end_date)
  f <- file.path(out_dir, paste0(out_prefix, "_layers.tif"))
  terra::writeRaster(r_layers, f, overwrite = TRUE)
  if (verbose) message("    已保存四层栅格: ", f)
  invisible(f)
}

## ============================================================================
##  S3 方法
## ============================================================================
print.rel_abundance <- function(x, ...) {
  cat("== eBird 相对丰度地图 ==\n")
  cat("  方法 :", x$method %||% "ebirdst", "\n")
  cat("  物种 :", x$species, "\n")
  cat("  时段 :", x$season, "\n")
  if (!is.null(x$r))
    cat("  栅格 :", paste(dim(x$r), collapse = " x "), "\n")
  if (!is.null(x$n_checklists))
    cat("  清单 :", x$n_checklists, "(建模", x$n_model %||% NA, ")\n")
  invisible(x)
}

plot.rel_abundance <- function(x, ...) {
  .plot_rel_abundance(x$r, x$region,
                      quantile_breaks = TRUE, n_quantiles = 10,
                      palette = "ebirdst",
                      map_title = paste0(x$species, " — ", x$season),
                      verbose = FALSE)
}


## ----------------------------------------------------------------------------
## 小样本去噪: 对连续层 (遭遇率/计数/丰度) 做轻度 focal 均值平滑。
## in_range 是 0/1 范围掩膜, 不平滑 (保持分布边界清晰); 只在研究区内部平滑。
## 数据量大、想完全复刻教程结果时可在 map_local(..., smooth = FALSE) 关闭。
## ----------------------------------------------------------------------------
.smooth_layers <- function(r, window = 3) {
  if (is.null(window) || !is.finite(window) || window <= 1) return(r)
  ## 注意: focal 配 mean 时窗口矩阵要用全 1 (若归一化, terra 会把 权重*值 再 mean, 等于多除一次窗口面积)
  w <- matrix(1, window, window)
  msk <- !is.na(r[[1]])
  smooth_one <- function(x) {
    s <- terra::focal(x, w, fun = function(z) mean(z, na.rm = TRUE))
    ## 边缘因窗口缺值产生 NA 的格子, 用原值回填, 保证研究区不出现空洞
    fill <- is.na(s) & !is.na(x)
    s[fill] <- x[fill]
    terra::mask(s, msk)
  }
  r[["encounter_rate"]] <- smooth_one(r[["encounter_rate"]])
  r[["count"]]          <- smooth_one(r[["count"]])
  r[["abundance"]]      <- r[["in_range"]] * r[["encounter_rate"]] *
                           r[["count"]]
  r
}
