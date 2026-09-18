## ============================================================================
##  ex.r —— f.rel.ab() 使用示例 (重点演示 method = "custom" 从零建模)
## ----------------------------------------------------------------------------
##  这个文件里的每一段都可以直接运行。建议在 RStudio 里打开本工程的 .Rproj,
##  这样工作目录就是工程根目录, 下面的相对路径 ("data", "output" 等) 才会生效。
##
##  【运行前只需准备一件事: 放数据】
##  把 eBird 基本数据集 (Custom Download, 勾选 "Include sampling event data")
##  解压得到的两个纯文本 (.txt/.tsv) 放进工程的 data/ 目录 (放子文件夹也行,
##  函数会自动递归查找):
##    * 清单文件 SED —— 文件名里带 "sampling"
##    * 观测文件 EBD —— 文件名以 "ebd" 开头但不带 "sampling"
##
##  【第一次运行会慢一些】
##  custom 法需要从网上免费下载环境栅格 (高程 + 土地覆盖, 无需密钥),
##  下载后会缓存到 data-raw/env-cache/, 之后同一区域重复运行直接读缓存。
##  缺少的 R 包会自动从 CRAN 安装。
##
##  【最省心的组合(强烈建议第一次先用它)】
##    region = list(type = "auto")   —— 研究区自动按你的数据范围来定
##    season = "year_round"          —— 不挑月份, 数据有啥用啥
##  这样几乎不会因为"地区/月份对不上"而失败。
## ============================================================================

## ---- 0. 加载函数 ----------------------------------------------------------
## 只需 source 一次; 函数会自己加载需要的 R 包, 你不必手动 library。
source("main.R")


## ============================================================================
##  一、method = "custom": 用你自己的 EBD/SED 从零建模 (本文件重点)
## ============================================================================

## ---- 示例 1【最稳, 第一次请先跑这个】---------------------------------------
## 研究区自动跟随数据范围, 时段用全年; 物种换成你数据里确实出现的鸟。
## 如果你不确定数据里有哪些鸟, 故意填一个错的名字, 报错信息会列出
## "记录最多的物种", 照抄一个即可。
ex1 <- f.rel.ab(
  species = "Northern Cardinal",        # 鸟的英文名(common name)
  season  = "year_round",               # 全年, 最不容易因月份对不上而失败
  region  = list(type = "auto"),        # 研究区 = 数据覆盖范围
  method  = "custom"
)

## ---- 示例 2: 指定季节(要和你数据的月份对得上) ------------------------------
## 可选季节: "breeding"(繁殖) / "nonbreeding"(越冬) /
##           "prebreeding_migration"(春季北迁前) /
##           "postbreeding_migration"(秋季南迁后) / "year_round"(全年) /
##           "weekly"(指定某一周) / "custom"(自定义起止日期)
## 若该季节在你的数据里没有记录, 函数会自动放宽到 "year_round" 并提示。
ex2 <- f.rel.ab(
  species = "Northern Cardinal",
  season  = "prebreeding_migration",    # 本例样例数据是 3 月, 对应春季迁徙前
  region  = list(type = "auto"),
  method  = "custom"
)

## ---- 示例 3: 一年中的第 N 周 (week 1-52) ----------------------------------
ex3 <- f.rel.ab(
  species = "Northern Cardinal",
  season  = "weekly",
  week    = 12,                         # 第 12 周(约 3 月中下旬)
  region  = list(type = "auto"),
  method  = "custom"
)

## ---- 示例 4: 完全自定义日期区间 -------------------------------------------
ex4 <- f.rel.ab(
  species     = "Northern Cardinal",
  season      = "custom",
  start_date  = "2025-03-05",
  end_date    = "2025-03-20",
  region      = list(type = "auto"),
  method      = "custom"
)

## ---- 示例 5: 指定行政区域而不是自动范围 ------------------------------------
## 注意: 你的数据必须落在该区域内; 若区域内记录太少, 函数会自动放宽
##       回"按数据范围出图"并给出中文提示, 保证你总能拿到一张图。
##
## 5a. 美国某州
ex5a <- f.rel.ab(
  species = "Northern Cardinal",
  season  = "year_round",
  region  = list(type = "state", name = "Alabama", country_iso = "US"),
  method  = "custom"
)

## 5b. 整个国家 (name 可写国名, 也可直接写 ISO2/ISO3 代码, 如 "US")
ex5b <- f.rel.ab(
  species = "Northern Cardinal",
  season  = "year_round",
  region  = list(type = "country", name = "United States of America"),
  method  = "custom"
)

## 5c. 美国本土 48 州 (不含阿拉斯加/夏威夷)
ex5c <- f.rel.ab(
  species = "Northern Cardinal",
  season  = "year_round",
  region  = list(type = "conus"),
  method  = "custom"
)

## 5d. 自定义多边形: 直接传一个 sf 对象 (例如用经纬度边界框)
##     下面框住 Alabama 西部一小块 (xmin,ymin,xmax,ymax, WGS84)
my_bbox <- sf::st_bbox(c(xmin = -87.8, ymin = 33.0, xmax = -87.0, ymax = 33.5),
                       crs = sf::st_crs(4326))
ex5d <- f.rel.ab(
  species = "Northern Cardinal",
  season  = "year_round",
  region  = sf::st_as_sfc(my_bbox),     # 直接传 sf/sfc 几何
  method  = "custom"
)

## ---- 示例 6: 用学名而不是英文名 -------------------------------------------
## 函数同时接受 common name 和学名(二名法)。
ex6 <- f.rel.ab(
  species = "Cardinalis cardinalis",    # Northern Cardinal 的学名
  season  = "year_round",
  region  = list(type = "auto"),
  method  = "custom"
)

## ---- 示例 7: 控制出图精细度 / 建模细节 -------------------------------------
ex7 <- f.rel.ab(
  species        = "Northern Cardinal",
  season         = "prebreeding_migration",
  region         = list(type = "auto"),
  method         = "custom",
  grid_res       = 3000,    # 预测网格边长(米); 调小=更细但更噪/更慢, 调大=更平滑
  elev_z         = 7,       # 高程下载精度等级(elevatr 的 z, 越大越精细)
  lc_agg         = 90,      # 土地覆盖聚合分辨率(米), 调大更快更平滑
  sample_radius  = 1500,    # 在每个点/网格周围多少米半径内统计环境(米)
  n_trees        = 500,     # 随机森林树的数量, 数据多可加大(更稳但更慢)
  smooth         = TRUE     # 对最终丰度做轻度空间平滑, 小样本去噪(默认就开)
)

## ---- 示例 8: 完全复刻教程, 关闭一切"兜底美化" ------------------------------
## smooth=FALSE 关闭空间平滑; 数据量大、想要原汁原味的模型输出时使用。
ex8 <- f.rel.ab(
  species = "Northern Cardinal",
  season  = "prebreeding_migration",
  region  = list(type = "auto"),
  method  = "custom",
  smooth  = FALSE
)

## ---- 示例 9: 数据不在默认 data/ 目录时, 手动指定路径 -----------------------
## 也可以直接用 data_dir 指向任意文件夹(会递归找 EBD/SED);
## 或用 ebd_file / sed_file 精确指定两个文件。
ex9 <- f.rel.ab(
  species  = "Northern Cardinal",
  season   = "year_round",
  region   = list(type = "auto"),
  method   = "custom",
  data_dir = "D:/my_ebird_data",        # 改成你自己的数据目录
  out_dir  = "my_maps"                  # 地图和栅格输出到这个文件夹
)


## ============================================================================
##  二、method = "ebirdst": 直接下载官方成品栅格 (快, 但需要 eBird 密钥)
##  下面这个 "yebsap-example" 是官方免密钥示例物种, 适合验证环境是否正常。
## ============================================================================
ex_ebirdst <- f.rel.ab(
  species    = "yebsap-example",
  season     = "breeding",
  region     = list(type = "state", name = "Michigan", country_iso = "US"),
  method     = "ebirdst",
  metric     = "median",     # median/mean/upper/lower 置信度
  resolution = "3km"         # 3km 或 27km
)


## ============================================================================
##  三、返回结果里有什么 / 怎么单独再画一次
## ============================================================================
## f.rel.ab 默认会: ① 弹出地图窗口; ② 在 out_dir(默认 "output/") 里保存
##   * PDF 地图
##   * <前缀>_abundance.tif  —— 最终相对丰度栅格(可在 GIS 软件打开)
##   * <前缀>_layers.tif     —— 四个分层栅格:
##       in_range         是否在分布范围内(0/1)
##       encounter_rate   校准后的遭遇率(0-1)
##       count            标准努力量下的期望个体数
##       abundance        = in_range × encounter_rate × count (最终相对丰度)
print(ex1)                  # 打印摘要
plot(ex1)                   # 重新画一次
terra::plot(ex1$r_layers[["encounter_rate"]])  # 只看某一层
ex1$diagnostics             # 模型诊断(阈值、测试集误差等)
ex1$n_checklists            # 用到多少条合格清单
ex1$season_used             # 实际使用的季节(若被自动放宽, 会和你填的不同)


## ============================================================================
##  四、常见问题 (排错)
## ============================================================================
## 1) 报"物种 ... 检测记录不足, 无法稳定建模"并列出一串鸟名:
##    说明你填的鸟在数据里太少或名字拼错; 从列出的常见鸟里抄一个英文名即可。
## 2) 提示"该季节记录不足, 已放宽到 year_round":
##    数据的月份和你选的季节对不上; 已自动用全年数据出图, 属正常保护机制。
## 3) 提示"研究区内记录不足, 已改用数据实际覆盖范围":
##    你选的行政区域和数据位置不符; 已自动按数据范围出图。
## 4) 第一次运行卡在"在线下载环境栅格": 网络较慢, 请耐心等待; 成功后会缓存。
## 5) 图太"碎"/椒盐: 小样本正常现象, 可增大 grid_res(如 5000) 或保持 smooth=TRUE。
## 6) 想固定随机结果以便复现: 设置 seed(默认 1), 或传入不同整数。
## ============================================================================
