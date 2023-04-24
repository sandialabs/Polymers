pub fn inverse_newton_raphson(y: &f64, f: &dyn Fn(&f64) -> f64, fp: &dyn Fn(&f64) -> f64, guess: &f64, &rel_tol: &f64, max_iters: &u8) -> f64
{
    let mut x = *guess;
    let mut y_minus_f: f64;
    let mut iters = 0;
    let mut residual_rel = 1.0;
    while residual_rel > rel_tol || &iters < max_iters
    {
        y_minus_f = y - f(&x);
        x += y_minus_f/fp(&x);
        iters += 1;
        residual_rel = (y_minus_f/y).abs();
    }
    x
}

pub fn inverse_langevin(y: &f64) -> f64
{
    if y <= &1e-3
    {
        3.0*y
    }
    else
    {
        inverse_newton_raphson(y, &|x: &f64| 1.0/x.tanh() - 1.0/x, &|x: &f64| 1.0/x.powi(2) - 1.0/x.sinh().powi(2), &((2.14234*y.powi(3) - 4.22785*y.powi(2) + 3.0*y)/(1.0 - y)/(0.71716*y.powi(3) - 0.41103*y.powi(2) - 0.39165*y + 1.0)), &1e-2, &100)
    }
}

pub fn integrate_1d(f: &dyn Fn(&f64) -> f64, x_min: &f64, x_max: &f64, num_points: &u128) -> f64
{
    let dx = (x_max - x_min)/(*num_points as f64);
    (0..*num_points).collect::<Vec::<u128>>().iter().map(|index| f(&(x_min + (0.5 + *index as f64)*dx))).sum::<f64>()*dx
}

pub fn integrate_1d_grid(f: &dyn Fn(&f64) -> f64, grid: &[f64], dx: &f64) -> f64
{
    grid.iter().map(|x| f(x)).sum::<f64>()*dx
}

pub fn integrate_2d(f: &dyn Fn(&f64, &f64) -> f64, x_min: &f64, x_max: &f64, y_min: &f64, y_max: &f64, num_points: &u128) -> f64
{
    let dx = (x_max - x_min)/(*num_points as f64);
    let dy = (y_max - y_min)/(*num_points as f64);
    let grid_x = (0..*num_points).collect::<Vec<u128>>().iter().map(|index| x_min + (0.5 + *index as f64)*dx).collect::<Vec<f64>>();
    let grid_y = (0..*num_points).collect::<Vec<u128>>().iter().map(|index| y_min + (0.5 + *index as f64)*dy).collect::<Vec<f64>>();
    grid_x.iter().flat_map(|x| grid_y.iter().map(|y| f(x, y))).sum::<f64>()*dx*dy
}

pub fn integrate_2d_symmetric(f: &dyn Fn(&f64, &f64) -> f64, x_min: &f64, x_max: &f64, num_points: &u128) -> f64
{
    let dx = (x_max - x_min)/(*num_points as f64);
    let grid = (0..*num_points).collect::<Vec<u128>>().iter().map(|index| x_min + (0.5 + *index as f64)*dx).collect::<Vec<f64>>();
    grid.iter().flat_map(|x| grid.iter().map(|y| f(x, y))).sum::<f64>()*dx.powi(2)
}

pub fn integrate_2d_grid(f: &dyn Fn(&f64, &f64) -> f64, grid_x: &[f64], grid_y: &[f64], da: &f64) -> f64
{
    grid_x.iter().flat_map(|x| grid_y.iter().map(|y| f(x, y))).sum::<f64>()*da
}

pub fn integrate_2d_symmetric_grid(f: &dyn Fn(&f64, &f64) -> f64, grid: &[f64], da: &f64) -> f64
{
    grid.iter().flat_map(|x| grid.iter().map(|y| f(x, y))).sum::<f64>()*da
}

pub fn lambert_w(x: &f64) -> f64
{
    let amount_of_iterations = ((x.log10()/3.0).ceil() as u8).max(4_u8);
    let mut w = 0.75*(x + 1.0).ln();
    for _ in 0..amount_of_iterations
    {
        w -= (w*w.exp() - x)/(w.exp()*(w + 1.0) - (w + 2.0)*(w*w.exp() - x)/(2.0 * w + 2.0));
    }
    w
}

pub fn bessel_i(nu: &u8, x: &f64) -> f64
{
    if nu == &0
    {
        bessel_i0(x)
    }
    else if nu == &1
    {
        bessel_i1(x)
    }
    else
    {
        -1.0
    }
}

fn bessel_i0(x: &f64) -> f64
{
    if x < &7.75
    {
        let coefficients = vec![
            1.0,
            2.499_999_999_999_999e-1,
            2.777_777_777_777_822e-2,
            1.736_111_111_110_237e-3,
            6.944_444_444_533_525e-5,
            1.929_012_345_132_199e-6,
            3.936_759_911_025_107e-8,
            6.151_186_727_044_392e-10,
            7.594_070_020_589_734e-12,
            7.593_897_933_698_363e-14,
            6.277_677_736_362_926e-16,
            4.347_097_041_532_722e-18,
            2.634_177_426_901_091e-20,
            1.139_430_377_448_228e-22,
            9.079_269_200_856_248e-25
        ];
        let t = 0.25*x.powi(2);
        1.0 + t*coefficients.iter().enumerate().map(|(i, c)| c*t.powi(i.try_into().unwrap())).sum::<f64>()
    }
    else if x < &500.0
    {
        let coefficients = vec![
             3.989_422_804_014_25e-1,
             4.986_778_506_049_619e-2,
             2.805_062_339_283_126e-2,
             2.922_112_251_660_478e-2,
             4.442_072_994_936_595e-2,
             1.309_705_746_058_567e-1,
            -3.350_522_802_317_27,
             2.330_257_115_835_147e2,
            -1.133_663_506_971_723e4,
             4.240_576_743_178_673e5,
            -1.231_570_285_956_987e7,
             2.802_319_381_552_675e8,
            -5.018_839_997_137_779e9,
             7.080_292_430_151_091e10,
            -7.842_610_821_248_111e11,
             6.768_257_378_540_965e12,
            -4.490_348_496_961_38e13,
             2.241_552_399_669_589e14,
            -8.134_264_678_656_593e14,
             2.023_910_973_916_877e15,
            -3.086_757_152_953_708e15,
             2.175_875_438_638_19e15
        ];
        x.exp()/x.sqrt()*coefficients.iter().enumerate().map(|(i, c)| c/x.powi(i.try_into().unwrap())).sum::<f64>()
    }
    else
    {
        let coefficients = vec![
            3.989_422_804_014_329e-1,
            4.986_778_504_914_345e-2,
            2.805_063_089_165_061e-2,
            2.921_790_968_539_151e-2,
            4.533_712_087_625_794e-2
        ];
        let expf = (0.5*x).exp();
        (expf/x.sqrt()*coefficients.iter().enumerate().map(|(i, c)| c/x.powi(i.try_into().unwrap())).sum::<f64>())*expf
    }
}

fn bessel_i1(x: &f64) -> f64
{
    if x < &7.75
    {
        let coefficients = vec![
            8.333_333_333_333_333e-2,
            6.944_444_444_444_341e-3,
            3.472_222_222_225_921e-4,
            1.157_407_407_354_987e-5,
            2.755_731_926_254_79e-7,
            4.920_949_692_800_671e-9,
            6.834_657_311_305_621e-11,
            7.593_969_849_687_574e-13,
            6.904_822_652_741_917e-15,
            5.220_157_095_351_373e-17,
            3.410_720_494_727_771e-19,
            1.625_212_890_947_171e-21,
            1.332_898_928_162_29e-23
        ];
        let t = 0.25*x.powi(2);
        let more_coefficients = vec![
            1.0,
            0.5,
            coefficients.iter().enumerate().map(|(i, c)| c*t.powi(i.try_into().unwrap())).sum::<f64>()
        ];
        0.5*x*more_coefficients.iter().enumerate().map(|(i, c)| c*t.powi(i.try_into().unwrap())).sum::<f64>()
    }
    else if x < &500.0
    {
        let coefficients = vec![
             3.989_422_804_014_406e-1,
            -1.496_033_551_613_111e-1,
            -4.675_104_253_598_537e-2,
            -4.090_895_951_581_637e-2,
            -5.719_036_414_430_205e-2,
            -1.528_189_554_374_492e-1,
             3.458_284_470_977_172e0,
            -2.426_181_371_595_021e2,
             1.178_785_865_993_44e4,
            -4.404_655_582_443_487e5,
             1.277_677_779_341_446e7,
            -2.903_390_398_236_656e8,
             5.192_386_898_222_206e9,
            -7.313_784_438_967_834e10,
             8.087_824_484_994_859e11,
            -6.967_602_516_005_787e12,
             4.614_040_809_616_582e13,
            -2.298_849_639_457_172e14,
             8.325_554_073_334_618e14,
            -2.067_285_045_778_906e15,
             3.146_401_654_361_325e15,
            -2.213_318_202_179_221e15
        ];
        x.exp()/x.sqrt()*coefficients.iter().enumerate().map(|(i, c)| c/x.powi(i.try_into().unwrap())).sum::<f64>()
    }
    else
    {
        let coefficients = vec![
             3.989_422_804_014_314e-1,
            -1.496_033_551_467_584e-1,
            -4.675_105_322_571_775e-2,
            -4.090_421_597_376_992e-2,
            -5.843_630_344_778_927e-2
        ];
        let expf = (0.5*x).exp();
        (expf/x.sqrt()*coefficients.iter().enumerate().map(|(i, c)| c/x.powi(i.try_into().unwrap())).sum::<f64>())*expf
    }
}

pub fn erf(x: &f64) -> f64
{
    1.0 - erfc(x)
}

pub fn erfc(x: &f64) -> f64
{
    erfcx(x)/(x.powi(2)).exp()
}

pub fn erfcx(x: &f64) -> f64
{
    if x >= &0.0
    {
        if x >= &50.0
        {
            if x > &5e7
            {
                0.564_189_583_547_756_3/x
            }
            else
            {
                0.564_189_583_547_756_3*(x.powi(2)*(x.powi(2) + 4.5) + 2.0)/(x*(x.powi(2)*(x.powi(2) + 5.0) + 3.75))
            }
        }
        else
        {
            erfcx_helper(&(400.0/(4.0 + x)))
        }
    }
    else if x < &-26.7
    {
        f64::MAX
    }
    else if x < &-6.1
    {
        2.0*(x.powi(2)).exp()
    }
    else
    {
        2.0*(x.powi(2)).exp() - erfcx_helper(&(400.0/(4.0 - x)))
    }
}

pub fn erfcx_helper(z: &f64) -> f64
{
    let zi = *z as u8;
    match zi
    {
        0_u8 =>
        {
            let t = 2.0*z - 1.0;
            7.087_803_245_410_644e-4 + (7.123_409_104_702_63e-4 + (3.577_907_729_759_774_2e-6 + (1.740_314_396_258_793_8e-8 + (8.171_066_004_730_779e-11 + (3.688_502_236_043_496e-13 + 1.591_703_855_111_111_2e-15*t)*t)*t)*t)*t)*t
        }
        1_u8 =>
        {
            let t = 2.0*z - 3.0;
            2.147_914_320_828_514_3e-3 + (7.268_640_236_737_999e-4 + (3.684_317_543_093_899_4e-6 + (1.807_184_127_214_92e-8 + (8.549_644_929_604_033e-11 + (3.885_203_751_853_429e-13 + 1.686_847_357_688_888_9e-15*t)*t)*t)*t)*t)*t
        }
        2_u8 =>
        {
            let t = 2.0*z - 5.0;
            3.616_525_593_563_017_5e-3 + (7.418_209_232_355_551e-4 + (3.794_831_995_752_824e-6 + (1.877_162_702_179_308_7e-8 + (8.948_471_512_241_509e-11 + (4.093_585_851_777_244e-13 + 1.787_206_146_488_889e-15*t)*t)*t)*t)*t)*t
        }
        3_u8 =>
        {
            let t = 2.0*z - 7.0;
            5.115_498_386_003_198e-3 + (7.572_284_073_479_166e-4 + (3.909_642_572_673_57e-6 + (1.950_416_870_430_047e-8 + (9.368_750_306_317_9e-11 + (4.314_392_595_907_966_5e-13 + 1.893_992_643_555_555_6e-15*t)*t)*t)*t)*t)*t
        }
        4_u8 =>
        {
            let t = 2.0*z - 9.0;
            6.645_751_317_267_305e-3 + (7.731_040_605_444_745e-4 + (4.028_951_058_939_944e-6 + (2.027_123_323_828_838_2e-8 + (9.811_763_132_170_91e-11 + (4.548_420_740_601_775e-13 + 2.007_635_221_333_333e-15*t)*t)*t)*t)*t)*t
        }
        5_u8 =>
        {
            let t = 2.0*z - 11.0;
            8.208_238_997_024_121e-3 + (7.894_662_961_188_171e-4 + (4.152_970_155_262_265e-6 + (2.107_469_334_454_465_7e-8 + (1.027_887_410_858_731_8e-10 + (4.796_520_139_061_334e-13 + 2.128_590_741_333_333_5e-15*t)*t)*t)*t)*t)*t
        }
        6_u8 =>
        {
            let t = 2.0*z - 13.0;
            9.803_953_727_535_219e-3 + (8.063_344_010_834_284e-4 + (4.281_924_132_973_699e-6 + (2.191_653_434_690_716_8e-8 + (1.077_153_513_656_547_1e-10 + (5.059_597_262_369_282e-13 + 2.257_346_268_444_444_6e-15*t)*t)*t)*t)*t)*t
        }
        7_u8 =>
        {
            let t = 2.0*z - 15.0;
            1.143_392_729_829_030_2e-2 + (8.237_285_838_319_657e-4 + (4.416_049_531_176_544e-6 + (2.279_886_142_621_198_7e-8 + (1.129_129_174_587_924e-10 + (5.338_618_936_581_688e-13 + 2.394_420_954_666_666_6e-15*t)*t)*t)*t)*t)*t
        }
        8_u8 =>
        {
            let t = 2.0*z - 17.0;
            1.309_923_287_881_465_4e-2 + (8.416_700_246_790_696e-4 + (4.555_595_898_845_751e-6 + (2.372_390_735_721_417_4e-8 + (1.183_978_932_660_269_6e-10 + (5.634_616_306_755_024e-13 + 2.540_367_964_444_444_6e-15*t)*t)*t)*t)*t)*t
        }
        9_u8 =>
        {
            let t = 2.0*z - 19.0;
            1.480_098_701_558_753_6e-2 + (8.601_809_294_634_594e-4 + (4.700_826_584_881_687e-6 + (2.469_404_076_019_731_5e-8 + (1.241_877_976_875_229_8e-10 + (5.948_689_037_032_026e-13 + 2.695_776_456_888_889e-15*t)*t)*t)*t)*t)*t
        }
        10_u8 =>
        {
            let t = 2.0*z - 21.0;
            1.654_035_173_939_407e-2 + (8.792_845_864_124_146e-4 + (4.852_019_579_300_175e-6 + (2.571_177_490_088_171e-8 + (1.303_012_853_423_082_1e-10 + (6.282_009_758_687_478e-13 + 2.861_273_735_111_111_2e-15*t)*t)*t)*t)*t)*t
        }
        11_u8 =>
        {
            let t = 2.0*z - 23.0;
            1.831_853_678_984_239_3e-2 + (8.990_054_264_789_172e-4 + (5.009_468_408_955_337e-6 + (2.677_977_707_421_807e-8 + (1.367_582_218_630_461_6e-10 + (6.635_828_774_535_271e-13 + 3.037_527_388_444_444_3e-15*t)*t)*t)*t)*t)*t
        }
        12_u8 =>
        {
            let t = 2.0*z - 25.0;
            2.013_680_196_421_427_7e-2 + (9.193_690_873_767_368e-4 + (5.173_483_091_410_427_6e-6 + (2.790_087_860_971_043_3e-8 + (1.435_797_640_280_904e-10 + (7.011_479_031_104_373e-13 + 3.225_247_6e-15*t)*t)*t)*t)*t)*t
        }
        13_u8 =>
        {
            let t = 2.0*z - 27.0;
            2.199_645_959_828_274_2e-2 + (9.404_024_815_536_678e-4 + (5.344_391_150_804_117e-6 + (2.907_808_553_804_937_5e-8 + (1.507_884_450_032_973e-10 + (7.410_381_364_749_92e-13 + 3.425_189_232e-15*t)*t)*t)*t)*t)*t
        }
        14_u8 =>
        {
            let t = 2.0*z - 29.0;
            2.389_887_718_722_632e-2 + (9.621_338_683_590_018e-4 + (5.522_538_699_804_901_5e-6 + (3.031_458_996_104_768_7e-8 + (1.584_082_649_729_633_4e-10 + (7.834_050_047_241_445e-13 + 3.638_155_356_444_444_5e-15*t)*t)*t)*t)*t)*t
        }
        15_u8 =>
        {
            let t = 2.0*z - 31.0;
            2.584_548_015_529_851_8e-2 + (9.845_929_306_782_012e-4 + (5.708_291_592_005_185e-6 + (3.161_378_216_916_483e-8 + (1.664_647_874_552_963e-10 + (8.284_098_592_878_54e-13 + 3.864_997_576_888_889e-15*t)*t)*t)*t)*t)*t
        }
        16_u8 =>
        {
            let t = 2.0*z - 33.0;
            2.783_775_478_347_469_8e-2 + (1.007_810_856_325_689_2e-3 + (5.902_036_649_379_221_6e-6 + (3.297_926_355_324_652e-8 + (1.749_852_415_926_845_7e-10 + (8.762_245_912_484_253e-13 + 4.106_620_648_888_889e-15*t)*t)*t)*t)*t)*t
        }
        17_u8 =>
        {
            let t = 2.0*z - 35.0;
            2.987_725_130_489_930_8e-2 + (1.031_820_424_505_735e-3 + (6.104_182_969_716_206e-6 + (3.441_486_035_954_272e-8 + (1.839_986_307_293_409e-10 + (9.270_322_736_636_504e-13 + 4.363_984_405_333_333_5e-15*t)*t)*t)*t)*t)*t
        }
        18_u8 =>
        {
            let t = 2.0*z - 37.0;
            3.196_558_717_859_645e-2 + (1.056_656_097_671_657_4e-3 + (6.315_163_319_241_458e-6 + (3.592_463_833_952_192e-8 + (1.935_358_475_878_117_3e-10 + (9.810_278_385_988_926e-13 + 4.638_106_081_777_777_6e-15*t)*t)*t)*t)*t)*t
        }
        19_u8 =>
        {
            let t = 2.0*z - 39.0;
            3.410_445_055_258_834e-2 + (1.082_354_119_135_053_2e-3 + (6.535_435_615_955_393e-6 + (3.751_291_834_853_352_4e-8 + (2.036_297_963_581_788_3e-10 + (1.038_418_783_303_728_1e-12 + 4.930_062_526_222_222e-15*t)*t)*t)*t)*t)*t
        }
        20_u8 =>
        {
            let t = 2.0*z - 41.0;
            3.629_560_392_829_242_5e-2 + (1.108_952_616_799_526_9e-3 + (6.765_484_509_551_836e-6 + (3.918_429_294_991_359e-8 + (2.143_155_220_213_377_5e-10 + (1.099_425_910_664_673_2e-12 + 5.240_994_910_222_222e-15*t)*t)*t)*t)*t)*t
        }
        21_u8 =>
        {
            let t = 2.0*z - 43.0;
            3.854_088_803_884_051e-2 + (1.136_491_713_417_542e-3 + (7.005_823_064_124_631e-6 + (4.094_364_408_371_858_4e-8 + (2.256_303_472_369_288_3e-10 + (1.164_284_101_136_199_3e-12 + 5.572_109_287_111_111e-15*t)*t)*t)*t)*t)*t
        }
        22_u8 =>
        {
            let t = 2.0*z - 45.0;
            4.084_222_595_478_596e-2 + (1.165_013_643_794_567_5e-3 + (7.256_994_550_234_3e-6 + (4.279_616_186_185_504e-8 + (2.376_140_171_100_502_3e-10 + (1.233_243_117_238_155_7e-12 + 5.924_680_236_444_444e-15*t)*t)*t)*t)*t)*t
        }
        23_u8 =>
        {
            let t = 2.0*z - 47.0;
            4.320_162_743_154_022e-2 + (1.194_562_879_391_727_1e-3 + (7.519_574_353_284_92e-6 + (4.474_736_455_396_099e-8 + (2.503_088_521_647_295e-10 + (1.306_568_440_030_047_7e-12 + 6.300_053_285_333_334e-15*t)*t)*t)*t)*t)*t
        }
        24_u8 =>
        {
            let t = 2.0*z - 49.0;
            4.562_119_351_381_047e-2 + (1.225_186_260_806_753e-3 + (7.794_172_005_555_192e-6 + (4.680_311_983_095_446e-8 + (2.637_599_098_397_842_6e-10 + (1.384_542_137_097_712e-12 + 6.699_647_740_444_444_6e-15*t)*t)*t)*t)*t)*t
        }
        25_u8 =>
        {
            let t = 2.0*z - 51.0;
            4.810_312_141_329_986_5e-2 + (1.256_933_138_643_219_5e-3 + (8.081_433_349_636_768e-6 + (4.896_966_733_568_202e-8 + (2.780_151_548_190_574_6e-10 + (1.467_463_761_160_988_5e-12 + 7.124_958_935_111_111e-15*t)*t)*t)*t)*t)*t
        }
        26_u8 =>
        {
            let t = 2.0*z - 53.0;
            5.064_970_967_698_334e-2 + (1.289_855_523_309_905_5e-3 + (8.382_042_841_456_88e-6 + (5.125_364_265_255_184e-8 + (2.931_256_384_967_550_7e-10 + (1.555_651_278_281_482_7e-12 + 7.577_560_782_222_223e-15*t)*t)*t)*t)*t)*t
        }
        27_u8 =>
        {
            let t = 2.0*z - 55.0;
            5.326_336_366_438_886_4e-2 + (1.324_008_244_325_697_5e-3 + (8.696_726_001_500_767e-6 + (5.366_210_275_039_68e-8 + (3.091_456_878_663_48e-10 + (1.649_442_024_082_849_4e-12 + 8.059_107_964_444_444e-15*t)*t)*t)*t)*t)*t
        }
        28_u8 =>
        {
            let t = 2.0*z - 57.0;
            5.594_660_135_350_001e-2 + (1.359_449_119_740_819e-3 + (9.026_252_023_301_638e-6 + (5.620_255_297_505_669_6e-8 + (3.261_331_041_050_314e-10 + (1.749_193_686_224_636_8e-12 + 8.571_338_168_888_888e-15*t)*t)*t)*t)*t)*t
        }
        29_u8 =>
        {
            let t = 2.0*z - 59.0;
            5.870_205_949_615_408_4e-2 + (1.396_239_136_322_364_7e-3 + (9.371_436_548_731_279e-6 + (5.888_297_567_026_528_5e-8 + (3.441_493_711_059_175_6e-10 + (1.855_285_310_975_186e-12 + 9.116_073_671_111_11e-15*t)*t)*t)*t)*t)*t
        }
        30_u8 =>
        {
            let t = 2.0*z - 61.0;
            6.153_250_014_514_477_5e-2 + (1.434_442_641_191_201_4e-3 + (9.733_144_620_101_681e-6 + (6.171_186_050_734_718e-8 + (3.632_598_741_829_53e-10 + (1.968_118_331_013_451_7e-12 + 9.695_223_84e-15*t)*t)*t)*t)*t)*t
        }
        31_u8 =>
        {
            let t = 2.0*z - 63.0;
            6.444_081_757_665_329e-2 + (1.474_127_545_638_313_2e-3 + (1.011_229_381_957_643_8e-5 + (6.469_823_660_593_325e-8 + (3.835_341_291_530_366_5e-10 + (2.088_117_611_438_512e-12 + 1.031_078_448e-14*t)*t)*t)*t)*t)*t
        }
        32_u8 =>
        {
            let t = 2.0*z - 65.0;
            6.743_004_563_313_039e-2 + (1.515_365_541_891_654e-3 + (1.050_985_760_688_832_9e-5 + (6.785_170_652_936_334e-8 + (4.050_460_219_481_114e-10 + (2.215_732_511_054_253_6e-12 + 1.096_484_211_555_555_5e-14*t)*t)*t)*t)*t)*t
        }
        33_u8 =>
        {
            let t = 2.0*z - 67.0;
            7.050_336_551_333_886e-2 + (1.558_232_333_649_571e-3 + (1.092_686_886_686_523e-5 + (7.118_248_223_961_351e-8 + (4.278_740_589_015_338_6e-10 + (2.351_437_952_227_442e-12 + 1.165_957_175_111_111_1e-14*t)*t)*t)*t)*t)*t
        }
        34_u8 =>
        {
            let t = 2.0*z - 69.0;
            7.366_411_403_794_46e-2 + (1.602_807_881_243_882e-3 + (1.136_442_367_877_820_8e-5 + (7.470_142_309_742_318e-8 + (4.521_016_277_747_649e-10 + (2.495_735_500_408_857e-12 + 1.239_723_825_777_777_7e-14*t)*t)*t)*t)*t)*t
        }
        35_u8 =>
        {
            let t = 2.0*z - 71.0;
            7.691_579_242_081_956e-2 + (1.649_176_662_344_788_9e-3 + (1.182_368_532_004_130_1e-5 + (7.842_007_599_378_154e-8 + (4.778_172_695_691_648e-10 + (2.649_154_440_381_572_5e-12 + 1.318_019_646_222_222_2e-14*t)*t)*t)*t)*t)*t
        }
        36_u8 =>
        {
            let t = 2.0*z - 73.0;
            8.026_207_557_809_462e-2 + (1.697_427_949_170_950_4e-3 + (1.230_588_851_730_989_1e-5 + (8.235_071_769_897_904e-8 + (5.051_149_610_985_711e-10 + (2.812_252_849_762_69e-12 + 1.401_088_963_555_555_5e-14*t)*t)*t)*t)*t)*t
        }
        37_u8 =>
        {
            let t = 2.0*z - 75.0;
            8.370_682_200_898_036e-2 + (1.747_656_103_221_265_7e-3 + (1.281_234_395_854_076_4e-5 + (8.650_639_951_503_644e-8 + (5.340_944_082_386_946e-10 + (2.985_618_662_088_755_5e-12 + 1.489_185_159_111_111e-14*t)*t)*t)*t)*t)*t
        }
        38_u8 =>
        {
            let t = 2.0*z - 77.0;
            8.725_408_428_446_171e-2 + (1.799_960_888_600_196_2e-3 + (1.334_444_308_008_949_3e-5 + (9.090_099_431_642_9e-8 + (5.648_613_497_261_646e-10 + (3.169_870_708_003_396e-12 + 1.582_569_779_555_555_6e-14*t)*t)*t)*t)*t)*t
        }
        39_u8 =>
        {
            let t = 2.0*z - 79.0;
            9.090_812_018_217_274e-2 + (1.854_447_805_065_77e-3 + (1.390_366_314_342_612e-5 + (9.554_924_606_254_991e-8 + (5.975_278_712_524_205e-10 + (3.365_659_736_609_91e-12 + 1.681_513_061_333_333_4e-14*t)*t)*t)*t)*t)*t
        }
        40_u8 =>
        {
            let t = 2.0*z - 81.0;
            9.467_340_450_807_549e-2 + (1.911_228_441_988_730_4e-3 + (1.449_157_261_654_500_5e-5 + (1.004_668_218_633_361_4e-7 + (6.322_127_295_979_1e-10 + (3.573_669_397_558_913e-12 + 1.786_293_159_111_111e-14*t)*t)*t)*t)*t)*t
        }
        41_u8 =>
        {
            let t = 2.0*z - 83.0;
            9.855_464_164_800_445e-2 + (1.970_420_854_472_562_2e-3 + (1.510_983_687_562_544_5e-5 + (1.056_703_666_767_598_4e-7 + (6.690_416_864_001_935e-10 + (3.794_617_185_082_434e-12 + 1.897_195_904e-14*t)*t)*t)*t)*t)*t
        }
        42_u8 =>
        {
            let t = 2.0*z - 85.0;
            1.025_567_788_947_009e-1 + (2.032_149_962_947_285_7e-3 + (1.576_022_424_296_218e-5 + (1.111_775_607_135_350_7e-7 + (7.081_478_511_009_766e-10 + (4.029_255_327_663_256e-12 + 2.014_514_307_555_555_6e-14*t)*t)*t)*t)*t)*t
        }
        43_u8 =>
        {
            let t = 2.0*z - 87.0;
            1.066_850_205_986_509_4e-1 + (2.096_547_977_614_873e-3 + (1.644_461_237_762_498_2e-5 + (1.170_071_796_202_615_3e-7 + (7.496_720_325_093_842e-10 + (4.278_371_618_608_592_5e-12 + 2.138_547_936e-14*t)*t)*t)*t)*t)*t
        }
        44_u8 =>
        {
            let t = 2.0*z - 89.0;
            1.109_448_431_938_644_4e-1 + (2.163_754_849_190_817e-3 + (1.716_499_503_571_965_6e-5 + (1.231_791_575_073_593_8e-7 + (7.937_630_983_149_963e-10 + (4.542_790_176_310_636e-12 + 2.269_602_565_333_333_3e-14*t)*t)*t)*t)*t)*t
        }
        45_u8 =>
        {
            let t = 2.0*z - 91.0;
            1.153_420_111_526_880_5e-1 + (2.233_918_747_454_642e-3 + (1.792_348_921_750_422_6e-5 + (1.297_146_528_824_599_7e-7 + (8.405_783_418_038_907e-10 + (4.823_372_120_641_802_5e-12 + 2.407_989_006_222_222_2e-14*t)*t)*t)*t)*t)*t
        }
        46_u8 =>
        {
            let t = 2.0*z - 93.0;
            1.198_825_939_268_409_5e-1 + (2.307_196_569_191_869e-3 + (1.872_234_271_895_893_7e-5 + (1.366_361_175_433_795_8e-7 + (8.902_838_548_849_328e-10 + (5.121_016_156_922_585e-12 + 2.554_022_711_111_111e-14*t)*t)*t)*t)*t)*t
        }
        47_u8 =>
        {
            let t = 2.0*z - 95.0;
            1.245_729_839_350_981_2e-1 + (2.383_754_477_180_957_6e-3 + (1.956_394_210_571_161e-5 + (1.439_673_684_773_947e-7 + (9.430_549_064_645_925e-10 + (5.436_659_058_313_422e-12 + 2.708_022_592e-14*t)*t)*t)*t)*t)*t
        }
        48_u8 =>
        {
            let t = 2.0*z - 97.0;
            1.294_199_156_614_243_8e-1 + (2.463_768_471_950_886e-3 + (2.045_082_112_747_588e-5 + (1.517_336_628_052_390_6e-7 + (9.990_763_250_638_903e-10 + (5.771_276_031_135_163e-12 + 2.870_309_955_555_555_5e-14*t)*t)*t)*t)*t)*t
        }
        49_u8 =>
        {
            let t = 2.0*z - 99.0;
            1.344_304_859_308_869_7e-1 + (2.547_424_998_108_082_3e-3 + (2.138_566_959_136_291_6e-5 + (1.599_617_757_990_044_2e-7 + (1.058_542_884_457_513_3e-9 + (6.125_880_953_678_788e-12 + 3.041_208_014_222_222e-14*t)*t)*t)*t)*t)*t
        }
        50_u8 =>
        {
            let t = 2.0*z - 101.0;
            1.396_121_754_343_456_2e-1 + (2.634_921_587_105_176_2e-3 + (2.237_134_271_257_256_8e-5 + (1.686_800_819_929_682_3e-7 + (1.121_659_691_044_499_7e-9 + (6.501_526_475_309_089e-12 + 3.221_039_450_666_667e-14*t)*t)*t)*t)*t)*t
        }
        51_u8 =>
        {
            let t = 2.0*z - 103.0;
            1.449_728_715_767_38e-1 + (2.726_467_538_398_244e-3 + (2.341_087_096_105_095e-5 + (1.779_186_393_952_637_8e-7 + (1.188_642_571_433_095_8e-9 + (6.899_303_966_505_428_4e-12 + 3.410_126_622_222_222_5e-14*t)*t)*t)*t)*t)*t
        }
        52_u8 =>
        {
            let t = 2.0*z - 105.0;
            1.505_208_927_277_461_9e-1 + (2.822_284_641_013_623_7e-3 + (2.450_747_042_271_339_8e-5 + (1.877_092_767_962_613_7e-7 + (1.259_718_458_758_337e-9 + (7.320_343_304_922_983e-12 + 3.608_788_904_888_888_7e-14*t)*t)*t)*t)*t)*t
        }
        53_u8 =>
        {
            let t = 2.0*z - 107.0;
            1.562_650_139_577_461e-1 + (2.922_607_937_619_662_7e-3 + (2.566_455_369_376_845e-5 + (1.980_856_841_565_446_2e-7 + (1.335_125_775_981_555_7e-9 + (7.765_812_489_104_676e-12 + 3.817_342_003_555_556e-14*t)*t)*t)*t)*t)*t
        }
        54_u8 =>
        {
            let t = 2.0*z - 109.0;
            1.622_144_943_462_073_8e-1 + (3.027_686_533_272_647_7e-3 + (2.688_574_132_653_456_3e-5 + (2.090_835_060_434_638_3e-7 + (1.415_114_814_424_073e-9 + (8.236_917_066_597_432e-12 + 4.036_095_745_777_778e-14*t)*t)*t)*t)*t)*t
        }
        55_u8 =>
        {
            let t = 2.0*z - 111.0;
            1.683_791_059_541_213e-1 + (3.137_784_451_079_308_3e-3 + (2.817_487_384_491_117_3e-5 + (2.207_404_380_704_578_2e-7 + (1.499_948_105_599_609e-9 + (8.734_899_366_193_081e-12 + 4.265_352_897_777_778e-14*t)*t)*t)*t)*t)*t
        }
        56_u8 =>
        {
            let t = 2.0*z - 113.0;
            1.747_691_645_565_937e-1 + (3.253_181_537_090_306_6e-3 + (2.953_602_434_734_436_5e-5 + (2.330_963_262_776_707_4e-7 + (1.589_900_784_358_244_5e-9 + (9.261_037_523_542_736e-12 + 4.505_407_310_222_222_4e-14*t)*t)*t)*t)*t)*t
        }
        57_u8 =>
        {
            let t = 2.0*z - 115.0;
            1.813_955_622_364_370_2e-1 + (3.374_174_416_809_7e-3 + (3.097_351_171_470_95e-5 + (2.461_932_693_759_229e-7 + (1.685_260_941_226_775_1e-9 + (9.816_644_294_285_49e-12 + 4.756_541_809_777_778e-14*t)*t)*t)*t)*t)*t
        }
        58_u8 =>
        {
            let t = 2.0*z - 117.0;
            1.882_698_019_444_366_5e-1 + (3.501_077_505_774_031_6e-3 + (3.249_191_444_001_427e-5 + (2.600_757_237_588_632e-7 + (1.786_329_961_738_837_7e-9 + (1.040_306_563_834_387_8e-11 + 5.019_026_583_111_111e-14*t)*t)*t)*t)*t)*t
        }
        59_u8 =>
        {
            let t = 2.0*z - 119.0;
            1.954_040_341_369_396_9e-1 + (3.634_224_076_721_132_6e-3 + (3.409_608_509_620_090_6e-5 + (2.747_906_111_701_763_6e-7 + (1.893_422_850_479_003_3e-9 + (1.102_167_907_532_359_9e-11 + 5.293_117_173_333_333e-14*t)*t)*t)*t)*t)*t
        }
        60_u8 =>
        {
            let t = 2.0*z - 121.0;
            2.028_110_956_065_188_6e-1 + (3.773_967_385_932_36e-3 + (3.579_116_545_759_241e-5 + (2.903_874_288_941_617_4e-7 + (2.006_868_537_484_9e-9 + (1.167_389_179_957_838_1e-11 + 5.579_052_309_333_333_5e-14*t)*t)*t)*t)*t)*t
        }
        61_u8 =>
        {
            let t = 2.0*z - 123.0;
            2.105_045_506_266_933_5e-1 + (3.920_681_861_392_565e-3 + (3.758_260_228_968_010_5e-5 + (3.069_183_623_188_687_7e-7 + (2.127_010_164_576_367_6e-9 + (1.236_113_855_106_289_9e-11 + 5.877_052_016e-14*t)*t)*t)*t)*t)*t
        }
        62_u8 =>
        {
            let t = 2.0*z - 125.0;
            2.184_987_345_370_333_3e-1 + (4.074_764_355_468_959e-3 + (3.947_616_382_098_671e-5 + (3.244_383_997_013_992e-7 + (2.254_205_349_151_868e-9 + (1.308_487_923_529_085_9e-11 + 6.187_315_326_222_222e-14*t)*t)*t)*t)*t)*t
        }
        63_u8 =>
        {
            let t = 2.0*z - 127.0;
            2.268_087_999_004_323e-1 + (4.236_635_464_862_852e-3 + (4.147_795_690_965_689_6e-5 + (3.430_054_489_450_281e-7 + (2.388_826_422_926_406_7e-9 + (1.384_659_629_281_851_4e-11 + 6.510_018_375_111_112e-14*t)*t)*t)*t)*t)*t
        }
        64_u8 =>
        {
            let t = 2.0*z - 129.0;
            2.354_507_653_698_870_4e-1 + (4.406_740_920_636_517e-3 + (4.359_444_491_622_47e-5 + (3.626_804_561_776_041_5e-7 + (2.531_260_643_085_32e-9 + (1.464_779_181_283_790_2e-11 + 6.845_312_263_111_111e-14*t)*t)*t)*t)*t)*t
        }
        65_u8 =>
        {
            let t = 2.0*z - 131.0;
            2.444_415_674_077_743_4e-1 + (4.585_553_051_160_578e-3 + (4.583_246_629_268_308_6e-5 + (3.835_275_259_003_303e-7 + (2.681_910_373_305_560_2e-9 + (1.548_998_439_088_475_8e-11 + 7.193_320_636_444_444e-14*t)*t)*t)*t)*t)*t
        }
        66_u8 =>
        {
            let t = 2.0*z - 133.0;
            2.537_991_150_063_426_7e-1 + (4.773_572_320_865_003e-3 + (4.819_925_389_653_418_5e-5 + (4.056_140_424_556_473_3e-7 + (2.841_193_232_087_116_4e-9 + (1.637_470_573_645_832e-11 + 7.554_137_982_222_222e-14*t)*t)*t)*t)*t)*t
        }
        67_u8 =>
        {
            let t = 2.0*z - 135.0;
            2.635_423_475_639_361_3e-1 + (4.971_328_947_708_378_5e-3 + (5.070_245_503_693_037e-5 + (4.290_107_925_426_818_5e-7 + (3.009_542_205_890_048e-9 + (1.730_349_702_534_734e-11 + 7.927_827_336_888_888e-14*t)*t)*t)*t)*t)*t
        }
        68_u8 =>
        {
            let t = 2.0*z - 137.0;
            2.736_912_960_773_234_5e-1 + (5.179_384_602_305_264e-3 + (5.335_015_225_832_66e-5 + (4.537_920_884_886_502e-7 + (3.187_405_724_581_438e-9 + (1.827_790_501_024_511e-11 + 8.314_418_236_444_444e-14*t)*t)*t)*t)*t)*t
        }
        69_u8 =>
        {
            let t = 2.0*z - 139.0;
            2.842_671_478_164_031_7e-1 + (5.398_334_191_669_514e-3 + (5.615_088_486_525_581e-5 + (4.800_358_919_649_474e-7 + (3.375_247_696_757_079_8e-9 + (1.929_947_788_808_346_8e-11 + 8.713_904_913_777_777e-14*t)*t)*t)*t)*t)*t
        }
        70_u8 =>
        {
            let t = 2.0*z - 141.0;
            2.952_923_146_534_852e-1 + (5.628_807_730_542_079_5e-3 + (5.911_367_118_991_331e-5 + (5.078_239_378_174_484e-7 + (3.573_547_502_585_171_4e-9 + (2.036_976_093_701_707e-11 + 9.126_244_261_333_333e-14*t)*t)*t)*t)*t)*t
        }
        71_u8 =>
        {
            let t = 2.0*z - 143.0;
            3.067_905_052_252_884e-1 + (5.871_472_303_274_54e-3 + (6.224_803_160_219_768e-5 + (5.372_418_576_620_094e-7 + (3.782_799_941_896_024e-9 + (2.149_029_193_044_454e-11 + 9.551_353_918_222_223e-14*t)*t)*t)*t)*t)*t
        }
        72_u8 =>
        {
            let t = 2.0*z - 145.0;
            3.187_868_011_117_332e-1 + (6.127_034_119_233_91e-3 + (6.556_401_225_970_764e-5 + (5.683_793_028_783_774e-7 + (4.003_515_135_339_238e-9 + (2.266_259_634_123_929_5e-11 + 9.989_110_976e-14*t)*t)*t)*t)*t)*t
        }
        73_u8 =>
        {
            let t = 2.0*z - 147.0;
            3.313_077_372_215_262_3e-1 + (6.396_240_664_679_808e-3 + (6.907_220_959_294_24e-5 + (6.013_300_666_188_594e-7 + (4.236_218_376_588_347e-9 + (2.388_818_234_707_369_7e-11 + 1.043_934_981_155_555_5e-13*t)*t)*t)*t)*t)*t
        }
        74_u8 =>
        {
            let t = 2.0*z - 149.0;
            3.443_813_865_804_133_4e-1 + (6.679_882_954_041_4e-3 + (7.278_379_551_860_356e-5 + (6.361_922_044_322_88e-7 + (4.481_449_933_651_445e-9 + (2.516_853_565_128_547_6e-11 + 1.090_186_138_311_111_1e-13*t)*t)*t)*t)*t)*t
        }
        75_u8 =>
        {
            let t = 2.0*z - 151.0;
            3.580_374_497_238_017_5e-1 + (6.978_797_883_488_269e-3 + (7.671_054_337_145_482e-5 + (6.730_681_530_891_739e-7 + (4.739_764_797_584_523e-9 + (2.650_511_414_114_305e-11 + 1.137_639_093_333_333_2e-13*t)*t)*t)*t)*t)*t
        }
        76_u8 =>
        {
            let t = 2.0*z - 153.0;
            3.723_073_489_011_972_7e-1 + (7.293_870_689_646_138e-3 + (8.086_485_454_267_072e-5 + (7.120_648_471_806_269e-7 + (5.011_732_376_974_588_4e-9 + (2.789_934_239_410_007_3e-11 + 1.186_263_761_422_222_2e-13*t)*t)*t)*t)*t)*t
        }
        77_u8 =>
        {
            let t = 2.0*z - 155.0;
            3.872_243_273_055_545e-1 + (7.626_037_516_254_98e-3 + (8.525_978_581_000_461e-5 + (7.532_938_330_517_133e-7 + (5.297_936_136_838_812e-9 + (2.935_260_605_416_409e-11 + 1.236_025_337_066_666_6e-13*t)*t)*t)*t)*t)*t
        }
        78_u8 =>
        {
            let t = 2.0*z - 157.0;
            4.028_235_535_461_694e-1 + (7.976_288_091_502_973e-3 + (8.990_907_734_243_825e-5 + (7.968_713_796_195_619e-7 + (5.598_973_180_736_040_5e-9 + (3.086_624_610_146_486_6e-11 + 1.286_884_194_666_666_8e-13*t)*t)*t)*t)*t)*t
        }
        79_u8 =>
        {
            let t = 2.0*z - 159.0;
            4.191_422_315_891_379e-1 + (8.345_668_518_695_046e-3 + (9.482_718_135_925_016e-5 + (8.429_185_856_178_314e-7 + (5.915_453_775_108_349e-9 + (3.244_155_303_434_747e-11 + 1.338_795_794_311_111e-13*t)*t)*t)*t)*t)*t
        }
        80_u8 =>
        {
            let t = 2.0*z - 161.0;
            4.362_197_163_946_378_6e-1 + (8.735_284_182_828_95e-3 + (1.000_292_914_206_68e-4 + (8.915_614_828_021_988e-7 + (6.248_000_815_078_86e-9 + (3.407_976_098_345_888e-11 + 1.391_710_717_688_888_8e-13*t)*t)*t)*t)*t)*t
        }
        81_u8 =>
        {
            let t = 2.0*z - 163.0;
            4.540_976_354_853_433e-1 + (9.146_302_775_554_824e-3 + (1.055_313_723_244_616_7e-4 + (9.429_311_346_463_863e-7 + (6.597_249_231_221_996e-9 + (3.578_204_179_547_656_4e-11 + 1.445_574_587_2e-13*t)*t)*t)*t)*t)*t
        }
        82_u8 =>
        {
            let t = 2.0*z - 165.0;
            4.728_200_166_851_233e-1 + (9.579_957_440_886_046e-3 + (1.113_501_905_800_006_7e-4 + (9.971_637_300_550_903e-7 + (6.963_845_336_995_697e-9 + (3.754_949_908_816_134_6e-11 + 1.500_328_071_288_889e-13*t)*t)*t)*t)*t)*t
        }
        83_u8 =>
        {
            let t = 2.0*z - 167.0;
            4.924_334_222_717_984e-1 + (1.003_755_004_390_949_7e-2 + (1.175_033_454_284_523_5e-4 + (1.054_400_671_618_896_7e-6 + (7.348_446_116_824_222e-9 + (3.938_316_232_643_575e-11 + 1.555_906_911_822_222_3e-13*t)*t)*t)*t)*t)*t
        }
        84_u8 =>
        {
            let t = 2.0*z - 169.0;
            5.129_870_897_920_926e-1 + (1.052_045_456_461_242_7e-2 + (1.240_093_003_749_499_7e-4 + (1.114_788_657_937_126_5e-6 + (7.751_718_455_056_87e-9 + (4.128_398_093_187_262_5e-11 + 1.612_241_968e-13*t)*t)*t)*t)*t)*t
        }
        85_u8 =>
        {
            let t = 2.0*z - 171.0;
            5.345_330_797_910_137e-1 + (1.103_012_061_880_072_7e-2 + (1.308_874_151_957_227e-4 + (1.178_479_759_537_451_5e-6 + (8.174_338_306_304_482e-9 + (4.325_281_844_951_708_4e-11 + 1.669_259_264e-13*t)*t)*t)*t)*t)*t
        }
        86_u8 =>
        {
            let t = 2.0*z - 173.0;
            5.571_264_307_116_93e-1 + (1.156_807_710_792_973_6e-2 + (1.381_579_783_803_665_2e-4 + (1.245_631_487_926_090_5e-6 + (8.616_989_807_896_932e-9 + (4.529_044_681_153_965e-11 + 1.726_880_108_444_444_3e-13*t)*t)*t)*t)*t)*t
        }
        87_u8 =>
        {
            let t = 2.0*z - 175.0;
            5.808_253_212_251_933e-1 + (1.213_593_599_950_387_8e-2 + (1.458_422_399_666_584e-4 + (1.316_406_857_309_571e-6 + (9.080_364_335_510_602e-9 + (4.739_754_071_312_462e-11 + 1.785_021_160_888_889e-13*t)*t)*t)*t)*t)*t
        }
        88_u8 =>
        {
            let t = 2.0*z - 177.0;
            6.056_912_402_529_337e-1 + (1.273_539_623_952_555e-2 + (1.539_624_447_225_886_4e-4 + (1.390_974_438_538_281_7e-6 + (9.565_159_503_230_623e-9 + (4.957_467_212_766_904e-11 + 1.843_594_556_444_444_4e-13*t)*t)*t)*t)*t)*t
        }
        89_u8 =>
        {
            let t = 2.0*z - 179.0;
            6.317_891_649_471_572e-1 + (1.336_824_779_828_703_2e-2 + (1.625_418_656_276_207_6e-4 + (1.469_508_404_833_405_5e-6 + (1.007_207_810_960_415_2e-8 + (5.182_230_499_568_071e-11 + 1.902_508_142_222_222_3e-13*t)*t)*t)*t)*t)*t
        }
        90_u8 =>
        {
            let t = 2.0*z - 181.0;
            6.591_877_468_972_532e-1 + (1.403_637_585_060_199_2e-2 + (1.716_048_376_025_970_7e-4 + (1.552_188_568_872_318_8e-6 + (1.060_182_703_153_528e-8 + (5.414_079_010_583_752e-11 + 1.961_665_514_666_666_7e-13*t)*t)*t)*t)*t)*t
        }
        91_u8 =>
        {
            let t = 2.0*z - 183.0;
            6.879_595_068_317_443e-1 + (1.474_176_509_136_586_8e-2 + (1.811_767_914_352_043_3e-4 + (1.639_200_410_823_058_4e-6 + (1.115_511_606_801_804_3e-8 + (5.653_036_019_492_569e-11 + 2.020_966_366_222_222_2e-13*t)*t)*t)*t)*t)*t
        }
        92_u8 =>
        {
            let t = 2.0*z - 185.0;
            7.181_810_380_872_997e-1 + (1.548_650_418_711_711_2e-2 + (1.912_842_878_455_092_4e-4 + (1.730_735_096_935_997_5e-6 + (1.173_265_673_611_360_8e-8 + (5.899_112_528_756_384e-11 + 2.080_306_533_333_333_4e-13*t)*t)*t)*t)*t)*t
        }
        93_u8 =>
        {
            let t = 2.0*z - 187.0;
            7.499_332_191_172_625e-1 + (1.627_279_036_404_478_3e-2 + (2.019_550_516_337_791_2e-4 + (1.826_989_488_320_334_8e-6 + (1.233_516_102_163_022_5e-8 + (6.152_306_831_216_908e-11 + 2.139_578_343_111_111_2e-13*t)*t)*t)*t)*t)*t
        }
        94_u8 =>
        {
            let t = 2.0*z - 189.0;
            7.833_014_353_128_349e-1 + (1.710_293_413_265_243e-2 + (2.132_180_058_506_332_8e-4 + (1.928_166_139_554_391_2e-6 + (1.296_334_008_735_434_2e-8 + (6.412_604_099_806_635e-11 + 2.198_670_894_222_222_3e-13*t)*t)*t)*t)*t)*t
        }
        95_u8 =>
        {
            let t = 2.0*z - 191.0;
            8.183_758_104_102_381e-1 + (1.797_936_414_904_422_3e-2 + (2.251_033_059_275_313e-4 + (2.034_473_286_801_817_5e-6 + (1.361_790_294_183_995e-8 + (6.679_976_008_397_248e-11 + 2.257_470_126_222_222e-13*t)*t)*t)*t)*t)*t
        }
        96_u8 =>
        {
            let t = 2.0*z - 193.0;
            8.552_514_477_568_512e-1 + (1.890_463_221_254_756e-2 + (2.376_423_737_037_125_5e-4 + (2.146_124_825_130_639e-6 + (1.429_955_507_187_052_3e-8 + (6.954_380_386_469_418e-11 + 2.315_859_368_888_888_7e-13*t)*t)*t)*t)*t)*t
        }
        97_u8 =>
        {
            let t = 2.0*z - 195.0;
            8.940_286_817_084_994e-1 + (1.988_141_839_912_72e-2 + (2.508_679_312_839_599_4e-4 + (2.263_340_274_758_523_3e-6 + (1.500_899_704_211_653e-8 + (7.235_760_907_504_394e-11 + 2.373_719_473_777_777_7e-13*t)*t)*t)*t)*t)*t
        }
        98_u8 =>
        {
            let t = 2.0*z - 197.0;
            9.348_133_394_287_079e-1 + (2.091_253_632_978_037e-2 + (2.648_140_346_599_848e-4 + (2.386_344_735_975_492_4e-6 + (1.574_692_306_547_218_3e-8 + (7.524_046_814_172_015e-11 + 2.430_929_127_111_111_4e-13*t)*t)*t)*t)*t)*t
        }
        99_u8 =>
        {
            let t = 2.0*z - 199.0;
            9.777_170_133_588_503e-1 + (2.200_093_857_283_048e-2 + (2.795_161_070_268_238e-4 + (2.515_368_832_524_531_6e-6 + (1.651_401_954_782_282e-8 + (7.819_152_682_936_823e-11 + 2.487_365_235_555_555_7e-13*t)*t)*t)*t)*t)*t
        }
        100_u8..=u8::MAX =>
        {
            1.0
        }
    }
}
