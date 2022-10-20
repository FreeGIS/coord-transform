use std::f64::consts::PI;
// 定义常量
const LL2MC: [[f64; 10]; 6] = [
    [
        -0.0015702102444,
        111320.7020616939,
        1704480524535203.0,
        -10338987376042340.0,
        26112667856603880.0,
        -35149669176653700.0,
        26595700718403920.0,
        -10725012454188240.0,
        1800819912950474.0,
        82.5,
    ],
    [
        0.0008277824516172526,
        111320.7020463578,
        647795574.6671607,
        -4082003173.641316,
        10774905663.51142,
        -15171875531.51559,
        12053065338.62167,
        -5124939663.577472,
        913311935.9512032,
        67.5,
    ],
    [
        0.00337398766765,
        111320.7020202162,
        4481351.045890365,
        -23393751.19931662,
        79682215.47186455,
        -115964993.2797253,
        97236711.15602145,
        -43661946.33752821,
        8477230.501135234,
        52.5,
    ],
    [
        0.00220636496208,
        111320.7020209128,
        51751.86112841131,
        3796837.749470245,
        992013.7397791013,
        -1221952.21711287,
        1340652.697009075,
        -620943.6990984312,
        144416.9293806241,
        37.5,
    ],
    [
        -0.0003441963504368392,
        111320.7020576856,
        278.2353980772752,
        2485758.690035394,
        6070.750963243378,
        54821.18345352118,
        9540.606633304236,
        -2710.55326746645,
        1405.483844121726,
        22.5,
    ],
    [
        -0.0003218135878613132,
        111320.7020701615,
        0.00369383431289,
        823725.6402795718,
        0.46104986909093,
        2351.343141331292,
        1.58060784298199,
        8.77738589078284,
        0.37238884252424,
        7.45,
    ],
];

const LLBAND: [f64; 6] = [75.0, 60.0, 45.0, 30.0, 15.0, 0.0];

const MCBAND: [f64; 6] = [
    12890594.86,
    8362377.87,
    5591021.0,
    3481989.83,
    1678043.12,
    0.0,
];
const MC2LL: [[f64; 10]; 6] = [
    [
        1.410526172116255e-8,
        0.00000898305509648872,
        -1.9939833816331,
        200.9824383106796,
        -187.2403703815547,
        91.6087516669843,
        -23.38765649603339,
        2.57121317296198,
        -0.03801003308653,
        17337981.2,
    ],
    [
        -7.435856389565537e-9,
        0.000008983055097726239,
        -0.78625201886289,
        96.32687599759846,
        -1.85204757529826,
        -59.36935905485877,
        47.40033549296737,
        -16.50741931063887,
        2.28786674699375,
        10260144.86,
    ],
    [
        -3.030883460898826e-8,
        0.00000898305509983578,
        0.30071316287616,
        59.74293618442277,
        7.357984074871,
        -25.38371002664745,
        13.45380521110908,
        -3.29883767235584,
        0.32710905363475,
        6856817.37,
    ],
    [
        -1.981981304930552e-8,
        0.000008983055099779535,
        0.03278182852591,
        40.31678527705744,
        0.65659298677277,
        -4.44255534477492,
        0.85341911805263,
        0.12923347998204,
        -0.04625736007561,
        4482777.06,
    ],
    [
        3.09191371068437e-9,
        0.000008983055096812155,
        0.00006995724062,
        23.10934304144901,
        -0.00023663490511,
        -0.6321817810242,
        -0.00663494467273,
        0.03430082397953,
        -0.00466043876332,
        2555164.4,
    ],
    [
        2.890871144776878e-9,
        0.000008983055095805407,
        -3.068298e-8,
        7.47137025468032,
        -0.00000353937994,
        -0.02145144861037,
        -0.00001234426596,
        0.00010322952773,
        -0.00000323890364,
        826088.5,
    ],
];

/** `bd2gcj` 将百度坐标系 (BD-09) 转换至 火星坐标系 (GCJ-02)
```
    let (x, y) = coord_transform::bd2gcj((118.0, 32.0));
    assert_eq!(x, 117.99349542486605);
    assert_eq!(y, 31.994068010063465);
```
*/
pub fn bd2gcj(bd: (f64, f64)) -> (f64, f64) {
    // 解构元组
    let (x, y) = bd;
    let x_pi = PI * 3000.0 / 180.0;
    let gcj_x = x - 0.0065;
    let gcj_y = y - 0.006;
    let z = (gcj_x.powi(2) + gcj_y.powi(2)).sqrt() - 0.00002 * ((gcj_y * x_pi).sin());
    let theta = gcj_y.atan2(gcj_x) - (gcj_x * x_pi).cos() * 0.000003;
    // 返回新元组
    (theta.cos() * z, theta.sin() * z)
}

/** `gcj2bd` 将火星坐标系 (GCJ-02) 转换至 百度坐标系 (BD-09)
```
    let (x, y) = coord_transform::gcj2bd((118.0, 32.0));
    assert_eq!(x, 118.00653128313961);
    assert_eq!(y, 32.00581846664105);
```
*/
pub fn gcj2bd(gcj: (f64, f64)) -> (f64, f64) {
    // 解构元组
    let (x, y) = gcj;
    let x_pi = PI * 3000.0 / 180.0;
    let z = (x.powi(2) + y.powi(2)).sqrt() + (y * x_pi).sin() * 0.00002;
    let theta = y.atan2(x) + (x * x_pi).cos() * 0.000003;
    // 返回新元组
    (theta.cos() * z + 0.0065, theta.sin() * z + 0.006)
}
/** `wgs2gcj` 将WGS84无偏移坐标 转换至 火星坐标系 (GCJ-02)
```
    let (x, y) = coord_transform::wgs2gcj((118.0, 32.0));
    assert_eq!(x, 118.00543101383846);
    assert_eq!(y, 31.997964381055795);
```
*/
pub fn wgs2gcj(wgs: (f64, f64)) -> (f64, f64) {
    let (x, y) = wgs;
    // 坐标在国外，直接返回
    if x < 72.004 || x > 137.8347 || y < 0.8293 || y > 55.8271 {
        return wgs;
    }
    let a = 6378245.0;
    let ee = 0.00669342162296594323;
    // 国内坐标
    let delta_x = x - 105.0;
    let delta_y = y - 35.0;

    let mut d_lat = -100.0
        + 2.0 * delta_x
        + 3.0 * delta_y
        + 0.2 * delta_y.powi(2)
        + 0.1 * delta_x * delta_y
        + 0.2 * (delta_x.abs().sqrt())
        + (20.0 * ((6.0 * delta_x * PI).sin()) + 20.0 * ((2.0 * delta_x * PI).sin())) * 2.0 / 3.0
        + (20.0 * ((delta_y * PI).sin()) + 40.0 * ((delta_y / 3.0 * PI).sin())) * 2.0 / 3.0
        + (160.0 * ((delta_y / 12.0 * PI).sin()) + 320.0 * ((delta_y * PI / 30.0).sin())) * 2.0
            / 3.0;

    let mut d_lon = 300.0
        + delta_x
        + 2.0 * delta_y
        + 0.1 * delta_x.powi(2)
        + 0.1 * delta_x * delta_y
        + 0.1 * delta_x.abs().sqrt()
        + (20.0 * ((6.0 * delta_x * PI).sin()) + 20.0 * ((2.0 * delta_x * PI).sin())) * 2.0 / 3.0
        + (20.0 * ((delta_x * PI).sin()) + 40.0 * ((delta_x / 3.0 * PI).sin())) * 2.0 / 3.0
        + (150.0 * ((delta_x / 12.0 * PI).sin()) + 300.0 * ((delta_x / 30.0 * PI).sin())) * 2.0
            / 3.0;

    let rad_lat = y / 180.0 * PI;
    let mut magic = rad_lat.sin();
    magic = 1.0 - ee * magic * magic;
    let magic_sqrt = magic.sqrt();

    d_lon = (d_lon * 180.0) / (a / magic_sqrt * (rad_lat.cos()) * PI);
    d_lat = (d_lat * 180.0) / ((a * (1.0 - ee)) / (magic * magic_sqrt) * PI);
    (x + d_lon, y + d_lat)
}
/** `gcj2wgs` 将火星坐标系 (GCJ-02)  转换至 WGS84无偏移坐标
```
    let (x, y) = coord_transform::gcj2wgs((118.0, 32.0));
    assert_eq!(x, 117.99456898616154);
    assert_eq!(y, 32.002035618944205);
```
*/
pub fn gcj2wgs(gcj: (f64, f64)) -> (f64, f64) {
    let (gcj_x, gcj_y) = gcj;
    let (x, y) = wgs2gcj(gcj);

    let d_lon = x - gcj_x;
    let d_lat = y - gcj_y;
    (gcj_x - d_lon, gcj_y - d_lat)
}

/** `bd2wgs` 将百度经纬度坐标系 (BD-09)  转换至 WGS84无偏移坐标
```
    let (x, y) = coord_transform::bd2wgs((118.0, 32.0));
    assert_eq!(x, 117.98808485929571);
    assert_eq!(y, 31.996121973877745);
```
*/
pub fn bd2wgs(bd: (f64, f64)) -> (f64, f64) {
    //百度先转火星，火星转84
    gcj2wgs(bd2gcj(bd))
}
/** `wgs2bd` 将WGS84无偏移坐标 转换至 百度经纬度坐标系 (BD-09)
```
    let (x, y) = coord_transform::wgs2bd((118.0, 32.0));
    assert_eq!(x, 118.01198481069936);
    assert_eq!(y, 32.00370423982076);
```
*/
pub fn wgs2bd(wgs: (f64, f64)) -> (f64, f64) {
    //百度先转火星，火星转84
    gcj2bd(wgs2gcj(wgs))
}

fn get_loop(mut lon: f64, min_v: f64, max_v: f64) -> f64 {
    while lon > max_v {
        lon -= max_v - min_v;
    }
    while lon < min_v {
        lon += max_v - min_v;
    }
    lon
}
fn get_range(mut lat: f64, min_v: f64, max_v: f64) -> f64 {
    lat = lat.max(min_v);
    lat = lat.min(max_v);
    lat
}
/** `bd_wgs2mkt` 将百度经纬度坐标(BD-09) 转换至 百度墨卡托坐标
```
    let (x, y) = coord_transform::bd_wgs2mkt((118.0, 32.0));
    assert_eq!(x, 13135842.840674074);
    assert_eq!(y, 3740459.445338771);
```
*/
pub fn bd_wgs2mkt(bdwgs: (f64, f64)) -> (f64, f64) {
    let (mut x, mut y) = bdwgs;
    x = get_loop(x, -180.0, 180.0);
    y = get_range(y, -74.0, 74.0);
    let mut cf: Option<[f64; 10]> = None;
    let mut i: usize = 0;
    for item in LLBAND {
        if y >= item {
            cf = Some(LL2MC[i]);
            break;
        }
        i += 1;
    }
    // 不为空匹配
    if let None = cf {
        for item in LLBAND {
            if y <= -item {
                cf = Some(LL2MC[i]);
                break;
            }
            i -= 1;
        }
    }
    //  let number = if a > 0 { 1 } else { -1 }; 类似js里的三元选择符
    match cf {
        Some(_cf)=> {
            x = _cf[0] + _cf[1] * x.abs();
            let cc = y.abs() / _cf[9];
            y = _cf[2]
                + _cf[3] * cc
                + _cf[4] * cc * cc
                + _cf[5] * cc * cc * cc
                + _cf[6] * cc * cc * cc * cc
                + _cf[7] * cc * cc * cc * cc * cc
                + _cf[8] * cc * cc * cc * cc * cc * cc;
            return (x.abs(), y.abs());
        }
        _ => {
            panic!("error occured");
        }
    }
}
/** `bd_mkt2wgs` 将百度墨卡托坐标 转换至 百度经纬度坐标(BD-09) 
```
    let (x, y) = coord_transform::bd_mkt2wgs((13135842.840674074,3740459.445338771));
    assert_eq!(x, 117.99999999999993);
    assert_eq!(y, 32.00000001036519);
```
*/
pub fn bd_mkt2wgs(bdmkt: (f64, f64)) -> (f64, f64) {
    let (mut x, mut y) = bdmkt;
    x = x.abs();
    y = y.abs();
    let mut cf: Option<[f64; 10]> = None;
    let mut i: usize = 0;
    for item in MCBAND {
        if y >= item {
            cf = Some(MC2LL[i]);
            break;
        }
        i += 1;
    }
    match cf {
        Some(_cf) => {
            x = _cf[0] + _cf[1] * (x.abs());
            let cc = y.abs() / _cf[9];
            y = _cf[2]
                + _cf[3] * cc
                + _cf[4] * cc * cc
                + _cf[5] * cc * cc * cc
                + _cf[6] * cc * cc * cc * cc
                + _cf[7] * cc * cc * cc * cc * cc
                + _cf[8] * cc * cc * cc * cc * cc * cc;
            return (x.abs(), y.abs());
        }
        _ => {
            panic!("error occured");
        }
    }
}

/** `wgs2bdmkt` 将无偏移的经纬度坐标 转换至 百度墨卡托坐标
```
    let (x, y) = coord_transform::wgs2bdmkt((118.0,32.0));
    assert_eq!(x, 13137176.998214714);
    assert_eq!(y, 3740943.32810748);
```
*/
pub fn wgs2bdmkt(wgs:(f64, f64)) -> (f64, f64) {
    bd_wgs2mkt(wgs2bd(wgs))
}
/** `bdmkt2wgs` 将百度墨卡托坐标 转换至 无偏移的经纬度坐标
```
    let (x, y) = coord_transform::bdmkt2wgs((13137176.998214714,3740943.32810748));
    assert_eq!(x, 117.9999832373631);
    assert_eq!(y, 31.999983711526742);
```
*/
pub fn bdmkt2wgs(bdmkt:(f64, f64)) -> (f64, f64) {
    bd2wgs(bd_mkt2wgs(bdmkt))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bd2gcj_test() {
        let (x, y) = bd2gcj((118.0, 32.0));
        assert_eq!(x, 117.99349542486605);
        assert_eq!(y, 31.994068010063465);
    }
    #[test]
    fn gcj2bd_test() {
        let (x, y) = gcj2bd((118.0, 32.0));
        assert_eq!(x, 118.00653128313961);
        assert_eq!(y, 32.00581846664105);
    }

    #[test]
    fn wgs2gcj_test() {
        let (x, y) = wgs2gcj((118.0, 32.0));
        assert_eq!(x, 118.00543101383846);
        assert_eq!(y, 31.997964381055795);
    }

    #[test]
    fn gcj2wgs_test() {
        let (x, y) = gcj2wgs((118.0, 32.0));
        assert_eq!(x, 117.99456898616154);
        assert_eq!(y, 32.002035618944205);
    }
    #[test]
    fn bd2wgs_test() {
        let (x, y) = bd2wgs((118.0, 32.0));
        assert_eq!(x, 117.98808485929571);
        assert_eq!(y, 31.996121973877745);
    }

    #[test]
    fn wgs2bd_test() {
        let (x, y) = wgs2bd((118.0, 32.0));
        assert_eq!(x, 118.01198481069936);
        assert_eq!(y, 32.00370423982076);
    }

    #[test]
    fn bd_wgs2mkt_test() {
        let (x, y) = bd_wgs2mkt((118.0, 32.0));
        assert_eq!(x, 13135842.840674074);
        assert_eq!(y, 3740459.445338771);
    }

    #[test]
    fn bd_mkt2wgs_test() {
        let (x, y) = bd_mkt2wgs((13135842.840674074,3740459.445338771));
        assert_eq!(x, 117.99999999999993);
        assert_eq!(y, 32.00000001036519);
    }

    #[test]
    fn wgs2bdmkt_test() {
        let (x, y) = wgs2bdmkt((118.0,32.0));
        assert_eq!(x, 13137176.998214714);
        assert_eq!(y, 3740943.32810748);
    }

    #[test]
    fn bdmkt2wgs_test() {
        let (x, y) = bdmkt2wgs((13137176.998214714,3740943.32810748));
        assert_eq!(x, 117.9999832373631);
        assert_eq!(y, 31.999983711526742);
    }
}