/// gamma function for integer n > 0
pub fn gamma(n: i32) -> f64 {
    if n < 1 {
        panic!("gamma not implemented for n < 1 (given value: n = {})", n);
    } else if ((n - 1) as usize) < GAMMA.len() {
        GAMMA[(n - 1) as usize]
    } else {
        std::f64::INFINITY
    }
}

// Table[Gamma[n], {n,1,171}]
const GAMMA: [f64; 171] = [
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0,
    362880.0              , 3.6288000000000000e006, 3.9916800000000000e007,
    4.7900160000000000e008, 6.2270208000000000e009, 8.7178291200000000e010,
    1.3076743680000000e012, 2.0922789888000000e013, 3.5568742809600000e014,
    6.4023737057280000e015, 1.2164510040883200e017, 2.4329020081766400e018,
    5.1090942171709440e019, 1.1240007277776077e021, 2.5852016738884977e022,
    6.2044840173323944e023, 1.5511210043330986e025, 4.0329146112660564e026,
    1.0888869450418352e028, 3.0488834461171386e029, 8.8417619937397020e030,
    2.6525285981219106e032, 8.2228386541779228e033, 2.6313083693369353e035,
    8.6833176188118865e036, 2.9523279903960414e038, 1.0333147966386145e040,
    3.7199332678990122e041, 1.3763753091226345e043, 5.2302261746660111e044,
    2.0397882081197443e046, 8.1591528324789773e047, 3.3452526613163807e049,
    1.4050061177528799e051, 6.0415263063373836e052, 2.6582715747884488e054,
    1.1962222086548019e056, 5.5026221598120889e057, 2.5862324151116818e059,
    1.2413915592536073e061, 6.0828186403426756e062, 3.0414093201713378e064,
    1.5511187532873823e066, 8.0658175170943879e067, 4.2748832840600256e069,
    2.3084369733924138e071, 1.2696403353658276e073, 7.1099858780486345e074,
    4.0526919504877217e076, 2.3505613312828786e078, 1.3868311854568984e080,
    8.3209871127413901e081, 5.0758021387722480e083, 3.1469973260387938e085,
    1.9826083154044401e087, 1.2688693218588416e089, 8.2476505920824707e090,
    5.4434493907744306e092, 3.6471110918188685e094, 2.4800355424368306e096,
    1.7112245242814131e098, 1.1978571669969892e100, 8.5047858856786232e101,
    6.1234458376886087e103, 4.4701154615126843e105, 3.3078854415193864e107,
    2.4809140811395398e109, 1.8854947016660503e111, 1.4518309202828587e113,
    1.1324281178206298e115, 8.9461821307829753e116, 7.1569457046263802e118,
    5.7971260207473680e120, 4.7536433370128417e122, 3.9455239697206587e124,
    3.3142401345653533e126, 2.8171041143805503e128, 2.4227095383672732e130,
    2.1077572983795277e132, 1.8548264225739844e134, 1.6507955160908461e136,
    1.4857159644817615e138, 1.3520015276784030e140, 1.2438414054641307e142,
    1.1567725070816416e144, 1.0873661566567431e146, 1.0329978488239059e148,
    9.9167793487094969e149, 9.6192759682482120e151, 9.4268904488832477e153,
    9.3326215443944153e155, 9.3326215443944153e157, 9.4259477598383594e159,
    9.6144667150351266e161, 9.9029007164861804e163, 1.0299016745145628e166,
    1.0813967582402909e168, 1.1462805637347084e170, 1.2265202031961379e172,
    1.3246418194518290e174, 1.4438595832024936e176, 1.5882455415227429e178,
    1.7629525510902447e180, 1.974506857221074e0182, 2.2311927486598136e184,
    2.5435597334721876e186, 2.9250936934930157e188, 3.3931086844518982e190,
    3.9699371608087209e192, 4.6845258497542907e194, 5.5745857612076059e196,
    6.6895029134491271e198, 8.0942985252734437e200, 9.8750442008336014e202,
    1.2146304367025330e205, 1.5061417415111409e207, 1.8826771768889261e209,
    2.3721732428800469e211, 3.0126600184576595e213, 3.8562048236258042e215,
    4.9745042224772874e217, 6.4668554892204737e219, 8.4715806908788205e221,
    1.1182486511960043e224, 1.4872707060906857e226, 1.9929427461615189e228,
    2.6904727073180505e230, 3.6590428819525487e232, 5.0128887482749917e234,
    6.9177864726194885e236, 9.6157231969410890e238, 1.3462012475717525e241,
    1.8981437590761710e243, 2.6953641378881628e245, 3.8543707171800728e247,
    5.5502938327393048e249, 8.0479260574719919e251, 1.1749972043909108e254,
    1.7272458904546389e256, 2.5563239178728656e258, 3.8089226376305697e260,
    5.7133839564458546e262, 8.6272097742332404e264, 1.3113358856834525e267,
    2.0063439050956824e269, 3.0897696138473509e271, 4.7891429014633939e273,
    7.4710629262828944e275, 1.1729568794264144e278, 1.8532718694937348e280,
    2.9467022724950383e282, 4.7147236359920613e284, 7.5907050539472187e286,
    1.2296942187394494e289, 2.0044015765453026e291, 3.2872185855342962e293,
    5.4239106661315888e295, 9.0036917057784374e297, 1.5036165148649990e300,
    2.5260757449731984e302, 4.2690680090047053e304, 7.2574156153079990e306
];

#[test]
fn test_gamma() {
    assert!(gamma(  1) == 1.0);
    assert!(gamma(  2) == 1.0);
    assert!(gamma(  3) == 2.0);
    assert!(gamma(170) == 4.2690680090047053e304);
    assert!(gamma(171) == 7.257415615307999e306);
    assert!(gamma(172).is_infinite());
}