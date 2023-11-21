#include <vector>
#include <cmath>

namespace data {
std::vector<double> dts = {0.2, 0.4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
std::vector<double> noise = {0.0, 3.52628601804259e-05, 0.00014947984792039608, 0.0003557569362467228, 0.000574050853713681, 0.0008010687876231075, 0.0010345987177312421, 0.0012728928570003401, 0.0015175881884682883, 0.0017637976446862192, 0.0020135697277017385, 0.0022653972995865484, 0.0025231914938045235, 0.0027794165373541003, 0.0030437654176523193, 0.003305867779530166, 0.003566952786824616, 0.0038330286045202507, 0.004098806518020103, 0.004368505041353395, 0.004636537001494376, 0.004906963983106297, 0.005181064460851167, 0.005459814320513834, 0.005731228969461196, 0.006005205037384404, 0.0062807635798247845, 0.006568777635662101, 0.006839555636779939, 0.007121599243086365, 0.0073993857278410085, 0.007681406791029645, 0.007963474798634514, 0.00824667948802273, 0.008527999754248573, 0.008817622484612405, 0.009100198837708833, 0.009385410406948088, 0.009674076157553613, 0.009957082096568705, 0.010246839853935069, 0.010542269449282776, 0.010823615217330204, 0.011115360297510592, 0.011409254187063518, 0.011703344176302588, 0.011990591528908462, 0.012295292495822487, 0.012582289994421238, 0.012870036304465413, 0.013179595508797827, 0.013451581938427694, 0.01374401979457549, 0.014047850960675352, 0.01434777464778414, 0.01464097565120466, 0.014952673145736457, 0.01524202406417912, 0.015533779053348214, 0.015823752513063175, 0.01611791951977581, 0.016429761432229607, 0.01674370891876689, 0.017030118950149536, 0.017332521380719803, 0.017638771326264825, 0.017933068649441133, 0.018258388035021524, 0.018535597084253865, 0.018851868135987156, 0.01912863613523952, 0.01944815458721639, 0.019767823575843026, 0.020058055572637894, 0.02036193356433419, 0.020672339265855044, 0.02097223790623132, 0.021293303453690654, 0.021614287985200365, 0.02190469640241079, 0.022179236602488075, 0.02253252204109979, 0.022815933500104708, 0.02313143595848678, 0.02345081928953667, 0.023749582390073787, 0.024019690037043784, 0.024396308028018187, 0.024706031114513072, 0.025002881834703086, 0.025295122647099735, 0.02560225525162434, 0.02594297581442481, 0.02625806940254291, 0.026551937707643886, 0.026876722159548606, 0.027159664337721063, 0.02752561373968782, 0.027850462500828135, 0.02813667597959067, 0.028454045798335218, 0.02874329978311689};

std::vector<std::vector <double>> Coeffs = {{-4.3504549427401994e-05, 0.00018996983968931822, 5.179654981227088e-06}, {-3.95897271932592e-05, 0.0001839581706814487, 8.272135904288011e-06, -5.16099776841178e-07, 1.833689921068164e-08, -2.640646501511814e-10}, {-4.2635391992097054e-05, 0.00018369597738974308, 8.921955239684975e-06, -6.071581106599217e-07, 1.8677959423439063e-08}, {-4.102177409106471e-05, 0.00018120416982138162, 1.0180341731412213e-05, -8.794837739058252e-07, 4.4901365676289396e-08, -9.233165017899216e-10}, {-8.427884108824172e-05, 0.00021083973352777608, 3.087661326898999e-06, -7.637603175812051e-08, 1.0782189985837739e-09, -6.101280885295214e-12}, {-0.00014397941095098963, 0.00023229072649020145, 9.930595861392308e-07}, {-3.924126845032596e-05, 0.00018076635044797694, 9.686076778276348e-06, -6.828873122654836e-07, 2.11483155880561e-08}, {-0.00018956122047978502, 0.0002382883769649114, 9.378131143201747e-07, -4.833764870762012e-09}, {-0.00018464655178129612, 0.0002357051355785347, 1.1328407794305888e-06, -9.620330271345615e-09, 3.6985452021968875e-11}, {-4.6381877607100126e-05, 0.00018721461028083216, 7.832267657865819e-06, -4.7259209311885137e-07, 1.2978901928392306e-08}, {-7.477079450400994e-05, 0.00020500372055766214, 3.9225658358390575e-06, -1.1526659908412697e-07, 1.4643285320845616e-09}, {-7.470413322214087e-05, 0.0002132643445915214, 2.3557766962725423e-06, -2.811730762651151e-08}, {-0.00012315510966583523, 0.00022081570593494051, 2.164216396083696e-06, -3.603441310022497e-08, 2.2658381604118924e-10, 9.118666022891758e-13}, {-0.00033535402640711906, 0.0002544207188724594, 4.286334196898537e-07}, {-0.00016516494411569597, 0.00023054520528897004, 1.456509981993346e-06, -1.6761733984385172e-08, 8.82145414697925e-11}, {-0.00019660210192317014, 0.00023534977698128225, 1.2146546744357314e-06, -1.1869067447618977e-08, 5.388179224376302e-11}, {-8.596526359498646e-05, 0.00020902943724014145, 3.5189336526516915e-06, -1.0767650782933352e-07, 1.99349768586106e-09, -1.553359334562716e-11}, {-3.0352073233267524e-05, 0.00018507848513959942, 7.151417903819826e-06, -3.569770464915893e-07, 9.887710104835302e-09, -1.0915940557272131e-10}, {-0.00012335858317680075, 0.0002268793981389805, 1.512901508267898e-06, -1.587545937364468e-08, 7.164403746223188e-11}, {-0.00022468610145360592, 0.00023809485264399694, 1.127316526272887e-06, -1.0801833131306316e-08, 4.982517557138517e-11}, {-0.0003876933965154244, 0.00025657588704816495, 4.1354871938798313e-07}, {-0.00021768747391224738, 0.00024132468763877807, 7.454137089276199e-07}, {-7.8969477970411e-05, 0.00020630941208504911, 3.805976996895073e-06, -1.1292564240772541e-07, 1.5163691036087063e-09}, {-0.0004844187309467271, 0.0002634792870520225, 3.0516067778516495e-07}, {-0.0002410517538964478, 0.0002402113267088728, 1.0452059676308838e-06, -9.90541654590542e-09, 5.6435795739738985e-11, -1.193499961015574e-13}, {-6.06010687854092e-05, 0.0002010162847179785, 4.2603684522716064e-06, -1.265532731094259e-07, 1.618305020038874e-09}, {-0.00021532428011694985, 0.00023946749597539386, 9.483168734129982e-07, -5.331121980867071e-09}, {-0.00042904478205683765, 0.00025992307018577024, 3.547667129558617e-07}, {-0.0002764645410975614, 0.00024524906453986885, 7.813885245960872e-07, -3.827496176453225e-09}, {-5.479609606657502e-05, 0.0002030954530497354, 3.6203310557820037e-06, -8.397074056640188e-08, 8.048608619281337e-10}, {-0.00023182138376774165, 0.0002397964661142527, 1.0275382709704844e-06, -8.686248038586552e-09, 3.5064536901914495e-11}, {-0.00031033742926044224, 0.0002485027925732776, 6.816325774103569e-07, -2.844772434227189e-09}, {0.00034207630102553194, 0.00011501083465393397, 1.1368364955628696e-05, -4.241430886733946e-07, 8.230209895475754e-09, -6.406270865064242e-11}, {-0.0003738409432642465, 0.0002562877632551079, 4.0948422885467474e-07}, {-0.0002539875035043552, 0.0002439022595797066, 7.994167243667845e-07, -3.822011861659083e-09}, {-0.00042269602364593435, 0.00025711840159663103, 4.760995019557176e-07, -1.3250801968874374e-09}, {-0.0004890388903855993, 0.00026230427479234934, 3.3354139905560586e-07}, {-0.00011988522739448271, 0.00022034094824465212, 2.1938278119349065e-06, -3.77210072772709e-08, 2.9157381643144654e-10}, {-5.296204727609291e-05, 0.0002176393793040466, 2.0839957712488383e-06, -3.283768994460494e-08, 3.0011656926547857e-10, -1.121197873720767e-12}, {-0.00015406392275099743, 0.0002294938962908306, 1.4187944140127802e-06, -1.199814650864976e-08}, {-0.002865162357342275, 0.00045758698077482546, -5.987714179175294e-06, 1.0128404513202884e-07, -8.036157464400697e-10, 2.487952024224264e-12}, {-0.0005316910896892516, 0.00026520966254753656, 2.90983749949382e-07}, {-0.00041874772353612914, 0.0002598704467222466, 3.495981471766633e-07}, {-0.00011913887252285385, 0.0002309149439545922, 1.1468639817494127e-06, -6.573259997662824e-09}, {-0.00045117698328634055, 0.0002573521690477333, 5.1633714028899e-07, -2.2481267486255083e-09, 5.3592319860303395e-12}, {0.00037724570031768406, 0.00012623232963043563, 9.069429117881773e-06, -2.804250182760144e-07, 4.438601023711654e-09, -2.7499512861955768e-11}, {-0.00042155438014845545, 0.00025747936806083937, 4.614794986546596e-07, -1.1952983284351445e-09}, {-0.0002169159258853014, 0.00024533196972672356, 5.798418368377113e-07}, {0.0055438053013385485, -0.000317210766182353, 2.1148272551390177e-05, -3.316259509477528e-07, 1.9798954810482776e-09}, {-0.0003623495778308297, 0.0002516677278466284, 6.47877582929673e-07, -3.58598180761647e-09, 1.0587960514228913e-11}, {-0.0003305679073209028, 0.00025373696912682904, 4.485688785805857e-07}, {-0.0004290503105051184, 0.00026014017266450516, 3.494546498803185e-07}, {-0.0004440799172992487, 0.00026137874134127807, 3.2978404727812093e-07}, {-0.00046938478257737073, 0.0002604538708565604, 4.047266441814135e-07, -8.587286235546852e-10}, {-0.00021920666638313574, 0.00023933187588139068, 1.0126219246739657e-06, -8.046742525625764e-09, 2.989824327857313e-11}, {-0.000591332614425271, 0.0002676116385320571, 2.6792836746952756e-07}, {-0.0003481843230590323, 0.00025224047846760536, 5.650052545184844e-07, -1.6920446950208605e-09}, {-0.00015961311871485855, 0.00023194860057080237, 1.299657163551894e-06, -1.2418158915408593e-08, 5.263201592307603e-11}, {-0.003287878524241769, 0.00036673417460461166, -6.4087096056379e-07}, {-0.00021532450808141833, 0.00023946752734497817, 9.483158511315388e-07, -5.33111233492802e-09}, {-0.000480521538112174, 0.0002647076405250958, 2.817262870231383e-07}, {-5.397524288514125, 0.4249123432416371, -0.013351388991809128, 0.0002096858989363652, -1.644796566037051e-06, 5.155122365818748e-09}, {-0.0002709678647765841, 0.00025198647889120636, 4.416836331478034e-07}, {-0.0006342793815902574, 0.0002697685001900074, 2.4421141051285105e-07}, {0.0020327771924538544, 0.00011742782662275745, 3.4428485942411454e-06, -2.923078494072297e-08, 9.8152718065011e-11}, {0.012200643084897309, -0.001052659221356987, 5.3979367715161156e-05, -1.0775077697737536e-06, 1.0667340470372605e-08, -4.1720692405582323e-11}, {-0.0006323272073794108, 0.0002694387967580393, 2.489642318777452e-07}, {-0.00010123126293212596, 0.00022028335586793705, 1.969148880285757e-06, -2.6777580402543572e-08, 1.550349326804841e-10}, {-0.00010465041378711175, 0.00022396125979665783, 1.6457254910647703e-06, -1.8169541142539557e-08, 8.495375393114166e-11}, {-0.00018421632603781747, 0.00023323235771846625, 1.3295183684361458e-06, -1.4299875437656938e-08, 7.132068325714401e-11}, {-0.00033603015162426107, 0.00025117961218080033, 5.801963178744245e-07, -1.1302009729472611e-09, -1.1372654788124432e-11}, {-7.52827987052316e-05, 0.00021498327409496284, 2.3434762107233437e-06, -3.938978318628905e-08, 3.60932515287072e-10, -1.2873894737423752e-12}, {-0.0005154578813732422, 0.0002648103125391961, 2.9394014970404696e-07}, {-0.0006790767603101072, 0.0002709225503343569, 2.3740530714776606e-07}, {1.4771521491267985, -0.07674513856322249, 0.0015039011537626185, -1.3035165575041164e-05, 4.233419447599374e-08}, {0.008094902692253025, -0.0003513826643891837, 1.7617823820713203e-05, -2.3915005229731073e-07, 1.623452995526419e-09, -4.356861980952508e-12}, {0.0011039909349001304, 0.00016533739582716932, 2.505735643586816e-06, -2.079662670460622e-08, 6.784812032873608e-11}, {-0.004372812667029078, 0.00044030435987532933, -2.3430734563720807e-06, 1.3068092851201322e-08}, {-0.000422948995340431, 0.0002610246831232043, 3.2692661838458874e-07}, {-1.1330781198326978, 0.07697415029721982, -0.0020727079728741337, 2.7939304941544395e-05, -1.87793870858854e-07, 5.035726345475158e-10}, {-0.002145323992150212, 0.0003822279632758504, -2.8832289126207003e-06, 3.829244571449628e-08, -1.7347541478847635e-10}, {-0.0002046256703360609, 0.00023641099124795814, 1.1724237123793624e-06, -1.1207069826132301e-08, 5.026340896875073e-11}, {-0.0009319597252704821, 0.00027855139403788385, 1.8064796719252016e-07}, {-0.0007087624138637921, 0.00027198987172826497, 2.285889913962418e-07}, {-0.00029320136230207884, 0.00024460041581467006, 8.953658891458359e-07, -7.170979303190041e-09, 2.895746293578936e-11}, {0.0010155872371839167, 0.00017102881528133535, 2.395190789156572e-06, -2.021513679119817e-08, 6.92689474234581e-11}, {-0.000334278770237853, 0.0002508352805737329, 6.159554020239741e-07, -2.307885041739073e-09}, {0.19217291026273198, -0.009959285634120931, 0.00020326945350034568, -1.7866938448468852e-06, 5.882516034251111e-09}, {58.77075263205767, -2.693866878086057, 0.04629921866439104, -0.0003535145272814109, 1.0119079128933063e-06}, {-2.644660181976289, 0.1278334256119563, -0.0023060051127038535, 1.8517778368705162e-05, -5.571699369429215e-08}, {-0.0001556969376702872, 0.0002282125150726725, 1.6507834898044813e-06, -2.387006723420216e-08, 2.059320950548885e-10, -7.215064621795049e-13}, {-0.00030006037387721, 0.0002496013835021259, 6.146560459136809e-07, -2.1131834413349235e-09}, {-0.0006322538275540038, 0.00026846481810418595, 2.784908891030938e-07, -2.1395965360323741e-10}, {-0.0006818216774061164, 0.0002709931343648663, 2.3730190562373946e-07}, {-0.00010054682596637166, 0.0002196619632007275, 2.0769074995303012e-06, -3.29547331401761e-08, 2.9171419876313334e-10, -1.0139071638953786e-12}, {-0.000297676847943934, 0.00026281807857465487, 2.794814783230791e-07}, {-0.0004080039067198109, 0.0002562368298003006, 4.927605064629074e-07, -1.4252742517606952e-09}, {-0.0006785819920644663, 0.00027081610566565077, 2.3935171581240387e-07}, {-0.00019349167625061503, 0.00023471354682907592, 1.2824232874628206e-06, -1.4625638628244959e-08, 9.942578236842674e-11, -2.6122805119072793e-13}, {-0.0003234416574072023, 0.00024835837525804915, 7.477619480593069e-07, -4.859340641255431e-09, 1.6404189741712833e-11}, {-0.0005686484761447223, 0.000267189361267373, 1.9487831313938428e-07, 2.9287397560874578e-09, -3.645336324891582e-11, 1.388971335221249e-13}, {-0.0002689172986106384, 0.00024325701551006467, 9.123204167387175e-07, -7.0001507806269006e-09, 2.563470300981467e-11}};

double approx(double dt, unsigned int id) {

    std::vector<double> coeffs = Coeffs[id];

    double ref = noise[id];
    double fit = 0;

    for (size_t i = 0; i < coeffs.size(); i++) {

        fit += coeffs[i] * std::pow(dt, i);

    }

    return fit/ref;

}

double approx_linear(double dt, unsigned int id) {

    std::vector<double> coeffs = Coeffs[id];

    double Ddt = dt - dts[id];

    double ref = noise[id];
    double fit = 0;

    if (id < noise.size()) { fit = noise[id] + (noise[id]-noise[id+1])*Ddt/(dts[id]-dts[id+1]);}
    else {fit = noise[id] + (noise[id]-noise[id-1])*Ddt/(dts[id]-dts[id-1]);}

    return fit/ref;

}

}