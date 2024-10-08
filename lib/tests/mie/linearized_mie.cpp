#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("LinMie construction", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();

    auto size_param = Eigen::VectorXd::LinSpaced(10, 0.1, 1.0);
    auto refractive_index = std::complex<double>(1.5, 0.0);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);

    auto result = mie.calculate(size_param, refractive_index, cos_angles, true);
}

TEST_CASE("LinMie construction2", "[sasktran2][mie]") {
    // testing the updates to An Bn code to ensure no Nans appear
    auto mie = sasktran2::mie::LinearizedMie();
    Eigen::VectorXd size_param(1000);
    size_param << 1.23816970e-04, 6.52382898e-04, 1.60330810e-03,
        2.97679442e-03, 4.77285419e-03, 6.99147586e-03, 9.63263963e-03,
        1.26963203e-02, 1.61824880e-02, 2.00911087e-02, 2.44221439e-02,
        2.91755509e-02, 3.43512829e-02, 3.99492890e-02, 4.59695140e-02,
        5.24118984e-02, 5.92763788e-02, 6.65628876e-02, 7.42713529e-02,
        8.24016986e-02, 9.09538447e-02, 9.99277067e-02, 1.09323196e-01,
        1.19140221e-01, 1.29378683e-01, 1.40038483e-01, 1.51119515e-01,
        1.62621670e-01, 1.74544833e-01, 1.86888889e-01, 1.99653715e-01,
        2.12839185e-01, 2.26445170e-01, 2.40471534e-01, 2.54918141e-01,
        2.69784846e-01, 2.85071505e-01, 3.00777966e-01, 3.16904073e-01,
        3.33449669e-01, 3.50414590e-01, 3.67798669e-01, 3.85601735e-01,
        4.03823611e-01, 4.22464118e-01, 4.41523073e-01, 4.61000287e-01,
        4.80895569e-01, 5.01208722e-01, 5.21939546e-01, 5.43087836e-01,
        5.64653385e-01, 5.86635979e-01, 6.09035402e-01, 6.31851433e-01,
        6.55083847e-01, 6.78732414e-01, 7.02796902e-01, 7.27277074e-01,
        7.52172687e-01, 7.77483498e-01, 8.03209255e-01, 8.29349706e-01,
        8.55904592e-01, 8.82873652e-01, 9.10256620e-01, 9.38053226e-01,
        9.66263195e-01, 9.94886251e-01, 1.02392211e+00, 1.05337049e+00,
        1.08323109e+00, 1.11350362e+00, 1.14418779e+00, 1.17528329e+00,
        1.20678982e+00, 1.23870706e+00, 1.27103470e+00, 1.30377242e+00,
        1.33691990e+00, 1.37047680e+00, 1.40444281e+00, 1.43881759e+00,
        1.47360079e+00, 1.50879207e+00, 1.54439109e+00, 1.58039750e+00,
        1.61681094e+00, 1.65363105e+00, 1.69085746e+00, 1.72848982e+00,
        1.76652775e+00, 1.80497087e+00, 1.84381882e+00, 1.88307119e+00,
        1.92272761e+00, 1.96278768e+00, 2.00325102e+00, 2.04411721e+00,
        2.08538587e+00, 2.12705658e+00, 2.16912892e+00, 2.21160249e+00,
        2.25447687e+00, 2.29775164e+00, 2.34142636e+00, 2.38550061e+00,
        2.42997395e+00, 2.47484594e+00, 2.52011615e+00, 2.56578412e+00,
        2.61184941e+00, 2.65831155e+00, 2.70517010e+00, 2.75242459e+00,
        2.80007456e+00, 2.84811953e+00, 2.89655903e+00, 2.94539258e+00,
        2.99461971e+00, 3.04423992e+00, 3.09425273e+00, 3.14465764e+00,
        3.19545416e+00, 3.24664179e+00, 3.29822002e+00, 3.35018834e+00,
        3.40254625e+00, 3.45529322e+00, 3.50842873e+00, 3.56195226e+00,
        3.61586329e+00, 3.67016128e+00, 3.72484569e+00, 3.77991599e+00,
        3.83537164e+00, 3.89121208e+00, 3.94743676e+00, 4.00404514e+00,
        4.06103666e+00, 4.11841074e+00, 4.17616684e+00, 4.23430437e+00,
        4.29282276e+00, 4.35172144e+00, 4.41099982e+00, 4.47065733e+00,
        4.53069337e+00, 4.59110735e+00, 4.65189868e+00, 4.71306676e+00,
        4.77461097e+00, 4.83653073e+00, 4.89882541e+00, 4.96149440e+00,
        5.02453708e+00, 5.08795284e+00, 5.15174104e+00, 5.21590106e+00,
        5.28043227e+00, 5.34533403e+00, 5.41060569e+00, 5.47624663e+00,
        5.54225618e+00, 5.60863369e+00, 5.67537852e+00, 5.74249001e+00,
        5.80996749e+00, 5.87781030e+00, 5.94601776e+00, 6.01458921e+00,
        6.08352398e+00, 6.15282137e+00, 6.22248071e+00, 6.29250132e+00,
        6.36288249e+00, 6.43362354e+00, 6.50472377e+00, 6.57618248e+00,
        6.64799896e+00, 6.72017250e+00, 6.79270240e+00, 6.86558794e+00,
        6.93882839e+00, 7.01242305e+00, 7.08637118e+00, 7.16067205e+00,
        7.23532493e+00, 7.31032908e+00, 7.38568377e+00, 7.46138826e+00,
        7.53744178e+00, 7.61384361e+00, 7.69059297e+00, 7.76768912e+00,
        7.84513130e+00, 7.92291874e+00, 8.00105067e+00, 8.07952633e+00,
        8.15834494e+00, 8.23750572e+00, 8.31700790e+00, 8.39685068e+00,
        8.47703329e+00, 8.55755493e+00, 8.63841480e+00, 8.71961212e+00,
        8.80114608e+00, 8.88301587e+00, 8.96522069e+00, 9.04775973e+00,
        9.13063217e+00, 9.21383720e+00, 9.29737400e+00, 9.38124174e+00,
        9.46543960e+00, 9.54996674e+00, 9.63482233e+00, 9.72000554e+00,
        9.80551553e+00, 9.89135145e+00, 9.97751246e+00, 1.00639977e+01,
        1.01508063e+01, 1.02379375e+01, 1.03253903e+01, 1.04131640e+01,
        1.05012575e+01, 1.05896702e+01, 1.06784010e+01, 1.07674492e+01,
        1.08568138e+01, 1.09464940e+01, 1.10364889e+01, 1.11267976e+01,
        1.12174191e+01, 1.13083527e+01, 1.13995974e+01, 1.14911524e+01,
        1.15830166e+01, 1.16751892e+01, 1.17676694e+01, 1.18604561e+01,
        1.19535485e+01, 1.20469457e+01, 1.21406466e+01, 1.22346505e+01,
        1.23289564e+01, 1.24235633e+01, 1.25184704e+01, 1.26136766e+01,
        1.27091811e+01, 1.28049829e+01, 1.29010810e+01, 1.29974746e+01,
        1.30941626e+01, 1.31911441e+01, 1.32884182e+01, 1.33859839e+01,
        1.34838402e+01, 1.35819861e+01, 1.36804208e+01, 1.37791432e+01,
        1.38781523e+01, 1.39774472e+01, 1.40770270e+01, 1.41768905e+01,
        1.42770369e+01, 1.43774651e+01, 1.44781741e+01, 1.45791631e+01,
        1.46804309e+01, 1.47819765e+01, 1.48837991e+01, 1.49858975e+01,
        1.50882707e+01, 1.51909178e+01, 1.52938377e+01, 1.53970295e+01,
        1.55004920e+01, 1.56042244e+01, 1.57082255e+01, 1.58124943e+01,
        1.59170298e+01, 1.60218310e+01, 1.61268969e+01, 1.62322263e+01,
        1.63378183e+01, 1.64436719e+01, 1.65497859e+01, 1.66561594e+01,
        1.67627913e+01, 1.68696804e+01, 1.69768259e+01, 1.70842266e+01,
        1.71918815e+01, 1.72997895e+01, 1.74079495e+01, 1.75163605e+01,
        1.76250214e+01, 1.77339311e+01, 1.78430886e+01, 1.79524927e+01,
        1.80621425e+01, 1.81720368e+01, 1.82821746e+01, 1.83925547e+01,
        1.85031761e+01, 1.86140376e+01, 1.87251382e+01, 1.88364768e+01,
        1.89480523e+01, 1.90598636e+01, 1.91719096e+01, 1.92841892e+01,
        1.93967012e+01, 1.95094447e+01, 1.96224183e+01, 1.97356211e+01,
        1.98490520e+01, 1.99627097e+01, 2.00765932e+01, 2.01907014e+01,
        2.03050331e+01, 2.04195873e+01, 2.05343627e+01, 2.06493582e+01,
        2.07645728e+01, 2.08800053e+01, 2.09956545e+01, 2.11115193e+01,
        2.12275985e+01, 2.13438911e+01, 2.14603959e+01, 2.15771116e+01,
        2.16940373e+01, 2.18111716e+01, 2.19285135e+01, 2.20460618e+01,
        2.21638154e+01, 2.22817730e+01, 2.23999336e+01, 2.25182959e+01,
        2.26368588e+01, 2.27556211e+01, 2.28745817e+01, 2.29937393e+01,
        2.31130929e+01, 2.32326411e+01, 2.33523830e+01, 2.34723171e+01,
        2.35924425e+01, 2.37127579e+01, 2.38332620e+01, 2.39539538e+01,
        2.40748321e+01, 2.41958955e+01, 2.43171430e+01, 2.44385734e+01,
        2.45601854e+01, 2.46819779e+01, 2.48039496e+01, 2.49260994e+01,
        2.50484261e+01, 2.51709283e+01, 2.52936051e+01, 2.54164550e+01,
        2.55394769e+01, 2.56626697e+01, 2.57860320e+01, 2.59095627e+01,
        2.60332606e+01, 2.61571244e+01, 2.62811529e+01, 2.64053449e+01,
        2.65296991e+01, 2.66542144e+01, 2.67788895e+01, 2.69037232e+01,
        2.70287142e+01, 2.71538613e+01, 2.72791633e+01, 2.74046190e+01,
        2.75302271e+01, 2.76559864e+01, 2.77818955e+01, 2.79079534e+01,
        2.80341587e+01, 2.81605103e+01, 2.82870068e+01, 2.84136470e+01,
        2.85404296e+01, 2.86673535e+01, 2.87944173e+01, 2.89216199e+01,
        2.90489598e+01, 2.91764360e+01, 2.93040472e+01, 2.94317920e+01,
        2.95596692e+01, 2.96876776e+01, 2.98158159e+01, 2.99440828e+01,
        3.00724771e+01, 3.02009976e+01, 3.03296428e+01, 3.04584116e+01,
        3.05873028e+01, 3.07163149e+01, 3.08454468e+01, 3.09746972e+01,
        3.11040648e+01, 3.12335484e+01, 3.13631466e+01, 3.14928582e+01,
        3.16226818e+01, 3.17526163e+01, 3.18826604e+01, 3.20128127e+01,
        3.21430720e+01, 3.22734369e+01, 3.24039063e+01, 3.25344788e+01,
        3.26651531e+01, 3.27959280e+01, 3.29268021e+01, 3.30577742e+01,
        3.31888430e+01, 3.33200071e+01, 3.34512653e+01, 3.35826163e+01,
        3.37140588e+01, 3.38455915e+01, 3.39772132e+01, 3.41089224e+01,
        3.42407179e+01, 3.43725985e+01, 3.45045627e+01, 3.46366094e+01,
        3.47687371e+01, 3.49009447e+01, 3.50332307e+01, 3.51655940e+01,
        3.52980331e+01, 3.54305469e+01, 3.55631339e+01, 3.56957929e+01,
        3.58285225e+01, 3.59613215e+01, 3.60941885e+01, 3.62271223e+01,
        3.63601215e+01, 3.64931848e+01, 3.66263109e+01, 3.67594985e+01,
        3.68927462e+01, 3.70260529e+01, 3.71594170e+01, 3.72928374e+01,
        3.74263127e+01, 3.75598417e+01, 3.76934228e+01, 3.78270550e+01,
        3.79607368e+01, 3.80944670e+01, 3.82282441e+01, 3.83620670e+01,
        3.84959342e+01, 3.86298444e+01, 3.87637965e+01, 3.88977889e+01,
        3.90318204e+01, 3.91658897e+01, 3.92999954e+01, 3.94341362e+01,
        3.95683109e+01, 3.97025180e+01, 3.98367563e+01, 3.99710244e+01,
        4.01053210e+01, 4.02396449e+01, 4.03739945e+01, 4.05083687e+01,
        4.06427662e+01, 4.07771855e+01, 4.09116253e+01, 4.10460845e+01,
        4.11805615e+01, 4.13150551e+01, 4.14495639e+01, 4.15840867e+01,
        4.17186221e+01, 4.18531688e+01, 4.19877254e+01, 4.21222906e+01,
        4.22568632e+01, 4.23914417e+01, 4.25260248e+01, 4.26606113e+01,
        4.27951998e+01, 4.29297889e+01, 4.30643774e+01, 4.31989639e+01,
        4.33335470e+01, 4.34681255e+01, 4.36026981e+01, 4.37372633e+01,
        4.38718199e+01, 4.40063666e+01, 4.41409020e+01, 4.42754248e+01,
        4.44099336e+01, 4.45444272e+01, 4.46789043e+01, 4.48133634e+01,
        4.49478032e+01, 4.50822225e+01, 4.52166200e+01, 4.53509942e+01,
        4.54853438e+01, 4.56196677e+01, 4.57539643e+01, 4.58882324e+01,
        4.60224707e+01, 4.61566778e+01, 4.62908525e+01, 4.64249933e+01,
        4.65590991e+01, 4.66931683e+01, 4.68271998e+01, 4.69611923e+01,
        4.70951443e+01, 4.72290545e+01, 4.73629217e+01, 4.74967446e+01,
        4.76305217e+01, 4.77642519e+01, 4.78979337e+01, 4.80315659e+01,
        4.81651471e+01, 4.82986760e+01, 4.84321513e+01, 4.85655717e+01,
        4.86989359e+01, 4.88322425e+01, 4.89654902e+01, 4.90986778e+01,
        4.92318039e+01, 4.93648672e+01, 4.94978664e+01, 4.96308002e+01,
        4.97636672e+01, 4.98964662e+01, 5.00291959e+01, 5.01618548e+01,
        5.02944418e+01, 5.04269556e+01, 5.05593947e+01, 5.06917580e+01,
        5.08240440e+01, 5.09562516e+01, 5.10883793e+01, 5.12204260e+01,
        5.13523902e+01, 5.14842708e+01, 5.16160663e+01, 5.17477755e+01,
        5.18793972e+01, 5.20109299e+01, 5.21423724e+01, 5.22737234e+01,
        5.24049816e+01, 5.25361458e+01, 5.26672145e+01, 5.27981866e+01,
        5.29290607e+01, 5.30598356e+01, 5.31905099e+01, 5.33210824e+01,
        5.34515518e+01, 5.35819167e+01, 5.37121760e+01, 5.38423283e+01,
        5.39723724e+01, 5.41023069e+01, 5.42321306e+01, 5.43618421e+01,
        5.44914403e+01, 5.46209239e+01, 5.47502915e+01, 5.48795419e+01,
        5.50086738e+01, 5.51376859e+01, 5.52665771e+01, 5.53953459e+01,
        5.55239911e+01, 5.56525116e+01, 5.57809059e+01, 5.59091728e+01,
        5.60373111e+01, 5.61653195e+01, 5.62931967e+01, 5.64209415e+01,
        5.65485527e+01, 5.66760289e+01, 5.68033689e+01, 5.69305714e+01,
        5.70576352e+01, 5.71845591e+01, 5.73113418e+01, 5.74379819e+01,
        5.75644784e+01, 5.76908300e+01, 5.78170353e+01, 5.79430932e+01,
        5.80690024e+01, 5.81947616e+01, 5.83203697e+01, 5.84458254e+01,
        5.85711274e+01, 5.86962745e+01, 5.88212655e+01, 5.89460992e+01,
        5.90707743e+01, 5.91952896e+01, 5.93196439e+01, 5.94438358e+01,
        5.95678643e+01, 5.96917281e+01, 5.98154260e+01, 5.99389567e+01,
        6.00623190e+01, 6.01855118e+01, 6.03085337e+01, 6.04313837e+01,
        6.05540604e+01, 6.06765626e+01, 6.07988893e+01, 6.09210391e+01,
        6.10430108e+01, 6.11648033e+01, 6.12864153e+01, 6.14078457e+01,
        6.15290932e+01, 6.16501567e+01, 6.17710349e+01, 6.18917267e+01,
        6.20122308e+01, 6.21325462e+01, 6.22526716e+01, 6.23726057e+01,
        6.24923476e+01, 6.26118958e+01, 6.27312494e+01, 6.28504070e+01,
        6.29693676e+01, 6.30881299e+01, 6.32066928e+01, 6.33250551e+01,
        6.34432157e+01, 6.35611733e+01, 6.36789269e+01, 6.37964752e+01,
        6.39138171e+01, 6.40309514e+01, 6.41478771e+01, 6.42645928e+01,
        6.43810976e+01, 6.44973902e+01, 6.46134694e+01, 6.47293342e+01,
        6.48449834e+01, 6.49604159e+01, 6.50756305e+01, 6.51906260e+01,
        6.53054014e+01, 6.54199556e+01, 6.55342873e+01, 6.56483955e+01,
        6.57622790e+01, 6.58759368e+01, 6.59893676e+01, 6.61025704e+01,
        6.62155441e+01, 6.63282875e+01, 6.64407995e+01, 6.65530791e+01,
        6.66651251e+01, 6.67769364e+01, 6.68885119e+01, 6.69998505e+01,
        6.71109511e+01, 6.72218127e+01, 6.73324340e+01, 6.74428141e+01,
        6.75529519e+01, 6.76628462e+01, 6.77724960e+01, 6.78819001e+01,
        6.79910576e+01, 6.80999674e+01, 6.82086283e+01, 6.83170392e+01,
        6.84251993e+01, 6.85331072e+01, 6.86407621e+01, 6.87481628e+01,
        6.88553083e+01, 6.89621975e+01, 6.90688293e+01, 6.91752028e+01,
        6.92813168e+01, 6.93871704e+01, 6.94927624e+01, 6.95980918e+01,
        6.97031577e+01, 6.98079589e+01, 6.99124944e+01, 7.00167633e+01,
        7.01207643e+01, 7.02244967e+01, 7.03279592e+01, 7.04311510e+01,
        7.05340709e+01, 7.06367180e+01, 7.07390913e+01, 7.08411896e+01,
        7.09430122e+01, 7.10445578e+01, 7.11458256e+01, 7.12468146e+01,
        7.13475236e+01, 7.14479518e+01, 7.15480982e+01, 7.16479618e+01,
        7.17475415e+01, 7.18468364e+01, 7.19458455e+01, 7.20445679e+01,
        7.21430026e+01, 7.22411485e+01, 7.23390048e+01, 7.24365705e+01,
        7.25338446e+01, 7.26308261e+01, 7.27275141e+01, 7.28239077e+01,
        7.29200058e+01, 7.30158076e+01, 7.31113121e+01, 7.32065183e+01,
        7.33014254e+01, 7.33960323e+01, 7.34903382e+01, 7.35843421e+01,
        7.36780430e+01, 7.37714402e+01, 7.38645326e+01, 7.39573193e+01,
        7.40497995e+01, 7.41419721e+01, 7.42338363e+01, 7.43253913e+01,
        7.44166360e+01, 7.45075696e+01, 7.45981911e+01, 7.46884998e+01,
        7.47784947e+01, 7.48681749e+01, 7.49575395e+01, 7.50465877e+01,
        7.51353185e+01, 7.52237312e+01, 7.53118248e+01, 7.53995984e+01,
        7.54870512e+01, 7.55741824e+01, 7.56609910e+01, 7.57474762e+01,
        7.58336373e+01, 7.59194732e+01, 7.60049832e+01, 7.60901664e+01,
        7.61750220e+01, 7.62595491e+01, 7.63437470e+01, 7.64276147e+01,
        7.65111515e+01, 7.65943565e+01, 7.66772290e+01, 7.67597680e+01,
        7.68419728e+01, 7.69238426e+01, 7.70053766e+01, 7.70865739e+01,
        7.71674338e+01, 7.72479554e+01, 7.73281380e+01, 7.74079808e+01,
        7.74874830e+01, 7.75666438e+01, 7.76454624e+01, 7.77239380e+01,
        7.78020700e+01, 7.78798574e+01, 7.79572996e+01, 7.80343957e+01,
        7.81111451e+01, 7.81875469e+01, 7.82636005e+01, 7.83393049e+01,
        7.84146596e+01, 7.84896638e+01, 7.85643167e+01, 7.86386175e+01,
        7.87125657e+01, 7.87861603e+01, 7.88594008e+01, 7.89322863e+01,
        7.90048162e+01, 7.90769898e+01, 7.91488062e+01, 7.92202649e+01,
        7.92913652e+01, 7.93621062e+01, 7.94324874e+01, 7.95025080e+01,
        7.95721673e+01, 7.96414647e+01, 7.97103995e+01, 7.97789709e+01,
        7.98471784e+01, 7.99150212e+01, 7.99824987e+01, 8.00496102e+01,
        8.01163550e+01, 8.01827325e+01, 8.02487421e+01, 8.03143830e+01,
        8.03796547e+01, 8.04445564e+01, 8.05090876e+01, 8.05732477e+01,
        8.06370359e+01, 8.07004516e+01, 8.07634943e+01, 8.08261633e+01,
        8.08884580e+01, 8.09503777e+01, 8.10119220e+01, 8.10730900e+01,
        8.11338814e+01, 8.11942953e+01, 8.12543314e+01, 8.13139889e+01,
        8.13732673e+01, 8.14321659e+01, 8.14906843e+01, 8.15488219e+01,
        8.16065780e+01, 8.16639520e+01, 8.17209436e+01, 8.17775519e+01,
        8.18337766e+01, 8.18896171e+01, 8.19450727e+01, 8.20001430e+01,
        8.20548274e+01, 8.21091254e+01, 8.21630364e+01, 8.22165600e+01,
        8.22696955e+01, 8.23224425e+01, 8.23748004e+01, 8.24267687e+01,
        8.24783469e+01, 8.25295345e+01, 8.25803311e+01, 8.26307360e+01,
        8.26807488e+01, 8.27303690e+01, 8.27795961e+01, 8.28284297e+01,
        8.28768692e+01, 8.29249142e+01, 8.29725641e+01, 8.30198186e+01,
        8.30666772e+01, 8.31131393e+01, 8.31592046e+01, 8.32048726e+01,
        8.32501428e+01, 8.32950148e+01, 8.33394881e+01, 8.33835623e+01,
        8.34272371e+01, 8.34705118e+01, 8.35133862e+01, 8.35558598e+01,
        8.35979321e+01, 8.36396028e+01, 8.36808715e+01, 8.37217377e+01,
        8.37622010e+01, 8.38022611e+01, 8.38419175e+01, 8.38811699e+01,
        8.39200178e+01, 8.39584610e+01, 8.39964989e+01, 8.40341312e+01,
        8.40713577e+01, 8.41081778e+01, 8.41445912e+01, 8.41805976e+01,
        8.42161966e+01, 8.42513879e+01, 8.42861711e+01, 8.43205459e+01,
        8.43545119e+01, 8.43880688e+01, 8.44212163e+01, 8.44539540e+01,
        8.44862816e+01, 8.45181989e+01, 8.45497054e+01, 8.45808009e+01,
        8.46114851e+01, 8.46417576e+01, 8.46716182e+01, 8.47010666e+01,
        8.47301025e+01, 8.47587255e+01, 8.47869355e+01, 8.48147321e+01,
        8.48421151e+01, 8.48690841e+01, 8.48956390e+01, 8.49217795e+01,
        8.49475052e+01, 8.49728160e+01, 8.49977116e+01, 8.50221918e+01,
        8.50462563e+01, 8.50699049e+01, 8.50931373e+01, 8.51159533e+01,
        8.51383527e+01, 8.51603353e+01, 8.51819009e+01, 8.52030492e+01,
        8.52237800e+01, 8.52440931e+01, 8.52639884e+01, 8.52834656e+01,
        8.53025246e+01, 8.53211651e+01, 8.53393870e+01, 8.53571900e+01,
        8.53745741e+01, 8.53915390e+01, 8.54080846e+01, 8.54242107e+01,
        8.54399172e+01, 8.54552039e+01, 8.54700706e+01, 8.54845172e+01,
        8.54985435e+01, 8.55121495e+01, 8.55253350e+01, 8.55380998e+01,
        8.55504439e+01, 8.55623670e+01, 8.55738692e+01, 8.55849502e+01,
        8.55956100e+01, 8.56058485e+01, 8.56156655e+01, 8.56250610e+01,
        8.56340349e+01, 8.56425870e+01, 8.56507174e+01, 8.56584258e+01,
        8.56657123e+01, 8.56725768e+01, 8.56790192e+01, 8.56850394e+01,
        8.56906374e+01, 8.56958132e+01, 8.57005666e+01, 8.57048976e+01,
        8.57088062e+01, 8.57122924e+01, 8.57153561e+01, 8.57179972e+01,
        8.57202159e+01, 8.57220119e+01, 8.57233854e+01, 8.57243363e+01,
        8.57248649e+01;
    auto refractive_index = std::complex<double>(1.4642667, 0.0);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);

    auto result = mie.calculate(size_param, refractive_index, cos_angles, true);
    REQUIRE(!std::isnan(result.values.Qext.sum()));
    REQUIRE(!std::isnan(result.values.Qsca.sum()));
    REQUIRE(!std::isnan(std::abs(result.values.S1.sum())));
    REQUIRE(!std::isnan(std::abs(result.values.S2.sum())));
}
/*
TEST_CASE("LinMie Dn test", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param = Eigen::VectorXd::LinSpaced(1, 62, 62);
    auto refractive_index = std::complex<double>(1.28, -1.37);
    int nstop = 50;
    Eigen::MatrixXcd Dn_matrix;
    mie.Dn(Dn_matrix, refractive_index, size_param, nstop);

    double ans_real = fabs(Dn_matrix(9, 0).real() - 0.004087);
    double ans_imag = fabs(Dn_matrix(9, 0).imag() - 1.0002620);
    REQUIRE(ans_real < 0.00001);
    REQUIRE(ans_imag < 0.00001);
}

TEST_CASE("LinMie An Bn test", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param = Eigen::VectorXd::LinSpaced(1, 50, 50);
    auto refractive_index = std::complex<double>(4.0 / 3.0, 0);
    double x = size_param.maxCoeff(); // largest size parameter, use to do the
                                      // calculations for highest N
    int N = int(x + 4.05 * pow(x, 0.33333) + 2.0) + 1;
    Eigen::MatrixXcd An_matrix;
    Eigen::MatrixXcd Bn_matrix;
    mie.An_Bn(refractive_index, size_param, N, An_matrix, Bn_matrix);

    REQUIRE(fabs(An_matrix(0, 0).real() - 0.5311058892948411929) < 0.00000001);
    REQUIRE(fabs(An_matrix(0, 0).imag() - 0.4990314856310943073) < 0.00000001);
    REQUIRE(fabs(Bn_matrix(0, 0).real() - 0.7919244759352004773) < 0.00001);
    REQUIRE(fabs(Bn_matrix(0, 0).imag() - 0.4059311522289938238) < 0.00001);

    auto size_param_2 = Eigen::VectorXd::LinSpaced(1, 2, 2);
    refractive_index = std::complex<double>(1.5, -1);
    x = size_param_2.maxCoeff(); // largest size parameter, use to do the
                                 // calculations for highest N
    N = int(x + 4.05 * pow(x, 0.33333) + 2.0) + 1;
    mie.An_Bn(refractive_index, size_param_2, N, An_matrix, Bn_matrix);

    REQUIRE(fabs(An_matrix(0, 0).real() - 0.5465202033970914511) < 0.00000001);
    REQUIRE(fabs(An_matrix(0, 0).imag() - 0.1523738572575972279) < 0.00000001);
    REQUIRE(fabs(Bn_matrix(0, 0).real() - 0.3897147278879423235) < 0.00001);
    REQUIRE(fabs(Bn_matrix(0, 0).imag() + 0.2278960752564908264) < 0.00001);

    refractive_index = std::complex<double>(1.1, -25);
    x = size_param_2.maxCoeff(); // largest size parameter, use to do the
                                 // calculations for highest N
    N = int(x + 4.05 * pow(x, 0.33333) + 2.0) + 1;
    mie.An_Bn(refractive_index, size_param_2, N, An_matrix, Bn_matrix);

    REQUIRE(fabs(An_matrix(1, 0).real() - 0.324433578437) < 0.0001);
    REQUIRE(fabs(An_matrix(1, 0).imag() - 0.465627763266) < 0.0001);
    REQUIRE(fabs(Bn_matrix(1, 0).real() - 0.060464399088) < 0.0001);
    REQUIRE(fabs(Bn_matrix(1, 0).imag() + 0.236805417045) < 0.0001);
}
*/

TEST_CASE("LinMie multiple size_params", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();

    auto size_param_og = Eigen::VectorXd::LinSpaced(2, 0.0004317428463995357,
                                                    0.0022637105717340984);
    auto refractive_index_og = std::complex<double>(1, -0.1);
    Eigen::VectorXd angles_og = Eigen::VectorXd::LinSpaced(15, 0, 180.0);
    Eigen::VectorXd cos_angles_og = cos(angles_og.array() * EIGEN_PI / 180.0);

    auto result_og =
        mie.calculate(size_param_og, refractive_index_og, cos_angles_og, true);

    auto size_param = Eigen::VectorXd::LinSpaced(2, 0.101, 0.2);
    auto refractive_index = std::complex<double>(0.75, 0.0);
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);

    auto result = mie.calculate(size_param, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext.size() - 2) < 0.0001);
    REQUIRE(fabs(result.values.Qsca.size() - 2) < 0.0001);

    REQUIRE(fabs(result.values.S1.cols() - 7) < 0.0001);
    REQUIRE(fabs(result.values.S1.rows() - 2) < 0.0001);
    REQUIRE(fabs(result.values.S2.cols() - 7) < 0.0001);
    REQUIRE(fabs(result.values.S2.rows() - 2) < 0.0001);
}

TEST_CASE("LinMie Qext Qsca test non absorbing (miepython)",
          "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    double lambda0 = 0.6328;
    double radius = 0.525;

    auto size_param_lambda = Eigen::VectorXd::LinSpaced(
        1, 2 * EIGEN_PI * radius / lambda0, 2 * EIGEN_PI * radius / lambda0);
    auto refractive_index = std::complex<double>(1.55, 0);
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);

    auto result =
        mie.calculate(size_param_lambda, refractive_index, cos_angles, true);

    REQUIRE(fabs(result.values.Qext(0) - 3.10543) < 0.00001);
    REQUIRE(fabs(result.values.Qsca(0) - 3.10543) < 0.00001);

    // MIEV0 Test Case 5
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_2 = Eigen::VectorXd::LinSpaced(1, 0.099, 0.099);
    result = mie.calculate(size_param_2, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000007) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.654225E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.654225E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.653815E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 1.574051E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.432172E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.652694E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 9.087788E-09) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 8.261265E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() + 1.651163E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 9.797186E-23) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() - 2.938374E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() + 1.649631E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 9.087788E-09) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() - 8.250360E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() + 1.648510E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 1.574051E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() - 1.427725E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() + 1.648100E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() - 1.648100E-04) < 1e-6);

    // MIEV0 Test Case 6
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_3 = Eigen::VectorXd::LinSpaced(1, 0.101, 0.101);
    result = mie.calculate(size_param_3, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000008) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 0.000008) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 2.048754E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.756419E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 2.048754E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.756419E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 2.048754E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.755965E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 1.774273E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.520629E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 2.048753E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.754726E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 1.024377E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 8.771198E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 2.048751E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() + 1.753033E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 1.845057E-15) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() - 3.247270E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 2.048750E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() + 1.751341E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.024375E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() - 8.759147E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 2.048749E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() + 1.750103E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 1.774269E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() - 1.515715E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 2.048749E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() + 1.749650E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 2.048749E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() - 1.749650E-04) < 1e-6);

    // MIEV0 Test Case 7
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_10 = Eigen::VectorXd::LinSpaced(1, 10, 10);
    result = mie.calculate(size_param_10, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.232265) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.580662E+01) <
            1e-5); // not enough digits for 1e-5
    REQUIRE(fabs(result.values.S1(0).imag() + 9.758097E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 5.580662E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(0).imag() + 9.758097E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() + 7.672879E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 1.087317E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() + 1.092923E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() - 9.629667E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 3.587894E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.756177E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 3.427411E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 8.082691E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() + 1.785905E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() + 5.232828E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() + 5.148748E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 7.027288E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 1.537971E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() + 8.329374E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 6.908338E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() - 2.152693E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() + 4.140427E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.876851E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() - 5.247557E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.923391E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() + 1.078568E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() + 3.608807E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() - 1.078568E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() - 3.608807E-02) < 1e-6);

    // MIEV0 Test Case 8
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_1000 = Eigen::VectorXd::LinSpaced(1, 1000, 1000);
    result = mie.calculate(size_param_1000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.997908) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 4.994770E+05) < 1e-1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.336502E+04) < 1e-2);
    REQUIRE(fabs(result.values.S2(0).real() - 4.994770E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.336502E+04) < 1e-2);

    REQUIRE(fabs(result.values.S1(1).real() + 3.999296E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(1).imag() + 3.316361E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(1).real() + 3.946018E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.147791E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(2).real() + 5.209852E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(2).imag() + 5.776614E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(2).real() + 1.970767E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(2).imag() + 6.937470E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(3).real() + 1.600887E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(3).imag() - 1.348013E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() + 4.152365E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).imag() - 1.143000E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(4).real() - 8.431720E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(4).imag() + 1.209493E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(4).real() + 4.261732E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).imag() - 5.535055E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 7.556092E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(5).imag() + 8.134810E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).real() - 4.218303E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).imag() - 9.100831E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(6).real() - 1.705778E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() - 4.842510E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(6).real() + 1.705778E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() + 4.842510E+02) < 1e-4);

    // OLD MIEV0 Test Case 1
    refractive_index = std::complex<double>(1.5, 0);
    result = mie.calculate(size_param_10, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.8820) < 1e-4);

    // OLD MIEV0 Test Case 2
    auto size_param_100 = Eigen::VectorXd::LinSpaced(1, 100, 100);
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.0944) < 1e-4);

    // OLD MIEV0 Test Case 3
    result = mie.calculate(size_param_1000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.0139) < 1e-4);

    // OLD MIEV0 Test Case 4
    auto size_param_5000 = Eigen::VectorXd::LinSpaced(1, 5000, 5000);
    result = mie.calculate(size_param_5000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.0086) < 1e-4);

    // test non-dialectric
    refractive_index = std::complex<double>(1.55, -0.1);
    result =
        mie.calculate(size_param_lambda, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 2.86165188243) < 1e-7);
    REQUIRE(fabs(result.values.Qsca(0) - 1.66424911991) < 1e-7);
}

TEST_CASE("LinMie Qext Qsca test absorbing (miepython)", "[sasktran2][mie]") {
    // MIEV0 Test Case 9
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param_1 = Eigen::VectorXd::LinSpaced(1, 1, 1);
    auto refractive_index = std::complex<double>(1.33, -0.00001);
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);
    auto result =
        mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.093923) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 9.395198E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 2.348800E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 2.281705E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 2.348800E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 2.281705E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 2.341722E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 2.217102E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 2.034830E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 1.938171E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 2.322408E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 2.046815E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 1.181704E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 1.075976E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 2.296081E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 1.828349E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 2.729533E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() - 6.702879E-03) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 2.269820E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 1.625401E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.114466E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 7.646326E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 2.250635E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.486170E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 1.942300E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.271557E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 2.243622E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 1.437106E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 2.243622E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 1.437106E-01) < 1e-6);

    // MIEV0 Test Case 10
    auto size_param_100 = Eigen::VectorXd::LinSpaced(1, 100, 100);
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.096594) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.101321E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.253302E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.243188E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(0).real() - 5.253302E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.243188E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(1).real() + 5.534573E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(1).imag() + 2.971881E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() + 8.467204E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.999470E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(2).real() - 1.710488E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.520096E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).real() - 3.310764E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).imag() + 2.709787E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() + 3.655758E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 8.769860E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() + 6.550512E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 4.675370E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 2.414318E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 5.380874E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() - 6.039011E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.169971E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 1.222996E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 3.283917E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).real() + 9.653812E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() - 1.474455E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(6).real() + 5.659205E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() - 4.650974E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).real() - 5.659205E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() + 4.650974E+01) < 1e-5);

    // // MIEV0 Test Case 11
    auto size_param_10000 = Eigen::VectorXd::LinSpaced(1, 10000, 10000);
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.723857) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.004089E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.010222E+07) < 1e1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.535815E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).real() - 5.010222E+07) < 1e1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.535815E+05) < 1e-1);

    REQUIRE(fabs(result.values.S1(1).real() - 3.786814E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(1).imag() + 7.654293E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).real() - 5.074755E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).imag() + 7.515986E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(2).real() + 2.731172E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(2).imag() - 1.326633E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).real() + 3.076558E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).imag() + 1.775975E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(3).real() + 1.061003E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(3).imag() + 1.930155E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() - 2.430920E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).imag() - 8.409836E+01) < 1e-5);

    REQUIRE(
        fabs(result.values.S1(4).real() + 1.05813972E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(4).imag() - 2.29841186E+01) <
        1e-5); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(4).real() - 5.90649314E+01) <
        1e-5); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(4).imag() + 5.37028316E+02) <
        1e-4); // different values using miepython with higher number of terms

    REQUIRE(
        fabs(result.values.S1(5).real() + 2.74885513E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(5).imag() - 2.29818128E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(5).real() + 8.03619653E+01) <
        1e-5); // different values using miepython with higher number of terms
    // the following criteria have been adjusted as confirmed by Dan Zawada - 6
    // digits matching with rounding
    REQUIRE(
        fabs(result.values.S2(5).imag() + 4.93916823E+00) <
        1e-5); // different values using miepython with higher number of terms

    REQUIRE(
        fabs(result.values.S1(6).real() + 1.82116215E+02) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(6).imag() + 9.51909674E+02) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).real() - 1.82116215E+02) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).imag() - 9.51909674E+02) <
        1e-3); // different values using miepython with higher number of terms

    // MIEV0 Test Case 12
    refractive_index = std::complex<double>(1.5, -1);
    auto size_param_55 = Eigen::VectorXd::LinSpaced(1, 0.055, 0.055);
    result = mie.calculate(size_param_55, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000011) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 1.014910E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 7.675259E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 8.343879E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 7.675259E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 8.343879E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 7.674331E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 8.343495E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 6.646948E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 7.225169E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 7.671794E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 8.342445E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 3.838246E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 4.169695E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 7.668328E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 8.341012E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 3.132066E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 2.037399E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 7.664863E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 8.339578E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 3.830082E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 4.171317E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 7.662326E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 8.338529E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 6.634986E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 7.221887E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 7.661398E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 8.338145E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 7.661398E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 8.338145E-05) < 1e-6);

    // MIEV0 Test Case 13
    auto size_param_56 = Eigen::VectorXd::LinSpaced(1, 0.056, 0.056);
    result = mie.calculate(size_param_56, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000012) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 1.033467E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 8.102381E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 8.807251E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 8.102381E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 8.807251E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 8.101364E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 8.806830E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 7.016844E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 7.626381E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 8.098587E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 8.805682E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 4.051865E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 4.401169E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 8.094795E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 8.804113E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 3.427277E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 2.229631E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 8.091003E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 8.802545E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 4.042932E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 4.402945E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 8.088228E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 8.801396E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 7.003755E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 7.622790E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 8.087213E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 8.800976E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 8.087213E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 8.800976E-05) < 1e-6);

    // MIEV0 Test Case 14
    result = mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.6634538) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.336321E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.840802E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 1.905153E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 5.840802E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 1.905153E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 5.657020E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 1.871997E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 5.001610E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 1.456112E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 5.175251E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 1.784426E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 2.879639E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 4.105398E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 4.563396E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 1.671665E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 3.622847E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 6.182646E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 4.002117E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 1.566427E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.748750E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.229586E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 3.621572E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.493910E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 3.056823E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.438460E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 3.488438E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 1.468286E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 3.488438E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 1.468286E-01) < 1e-6);

    // MIEV0 Test Case 15
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.283697) < 1e-3);
    REQUIRE(fabs(result.values.Qext(0) - 2.097502) < 1e-2);

    REQUIRE(fabs(result.values.S1(0).real() - 5.243754E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(0).imag() + 2.934167E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(0).real() - 5.243754E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(0).imag() + 2.934167E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(1).real() - 4.049055E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.898456E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() - 2.019198E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() - 3.110731E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() + 2.646835E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.929564E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).real() - 9.152743E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 7.470202E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 1.268890E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(3).imag() - 2.397474E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).real() + 1.232914E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).imag() + 7.823167E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 5.149886E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 2.290736E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).real() + 7.173357E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.655464E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 1.605395E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.418642E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).real() - 1.448052E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.393594E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(6).real() + 2.029360E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() - 4.384435E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() - 2.029360E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() + 4.384435E+00) < 1e-6);

    // MIEV0 Test Case 16
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.236575) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.004368) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.010919E+07) < 1e1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.753404E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).real() - 5.010919E+07) < 1e1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.753404E+05) < 1e-1);

    REQUIRE(fabs(result.values.S1(1).real() + 3.690394E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.573897E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).real() + 9.333175E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.839736E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(2).real() - 2.391551E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(2).imag() - 3.247786E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).real() + 1.202951E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).imag() + 1.899647E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(3).real() + 2.607463E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(3).imag() - 7.414859E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() - 1.013073E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(3).imag() + 1.064666E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(4).real() + 6.183154E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(4).imag() - 2.264970E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(4).real() - 1.334826E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.800859E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(5).real() + 3.368019E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(5).imag() - 2.115750E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(5).real() - 2.293862E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.996754E+03) < 1e-3);

    REQUIRE(
        fabs(result.values.S1(6).real() + 2.18471714E+02) <
        1e-4); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(6).imag() + 2064.61012226) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).real() - 2.18471714E+02) <
        1e-4); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).imag() - 2064.61012226) <
        1e-3); // different values using miepython with higher number of terms

    // MIEV0 Test Case 17
    refractive_index = std::complex<double>(10, -10);
    result = mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.049405) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.532993E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 6.332483E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 4.179305E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 6.332483E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 4.179305E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 6.162264E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 4.597163E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 5.573186E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 2.954338E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 5.736317E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 5.602514E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 3.525107E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 5.921611E-03) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 5.238628E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 6.675352E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 7.881172E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 3.435544E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 4.825816E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 7.434033E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.881212E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 6.028739E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 4.570214E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 7.809867E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 3.793898E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 7.473279E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 4.485464E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 7.912365E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 4.485464E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 7.912365E-01) < 1e-6);

    // MIEV0 Test Case 18
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.836785) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.071124E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.177811E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(0).imag() + 2.633811E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(0).real() - 5.177811E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(0).imag() + 2.633811E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(1).real() - 5.227436E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.270012E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() + 2.380252E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() + 3.872567E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() + 2.705712E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(2).imag() + 3.951751E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).real() - 2.585821E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).imag() - 3.323624E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(3).real() - 1.008860E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 4.663027E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).real() + 3.479935E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 4.364245E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(4).real() + 1.505640E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(4).imag() - 4.333057E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).real() - 1.360634E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).imag() + 4.238302E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 4.510770E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(5).imag() - 5.199554E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() - 4.474564E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).imag() + 5.452513E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() + 4.145383E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() + 1.821808E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).real() - 4.145383E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() - 1.821808E+01) < 1e-5);

    // MIEV0 Test Case 19
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.795393) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.005914E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.014786E+07) < 1e1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.206004E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).real() - 5.014786E+07) < 1e1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.206004E+05) < 1e-1);

    REQUIRE(fabs(result.values.S1(1).real() + 4.080090E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(1).imag() + 2.664399E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).real() - 3.351286E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).imag() - 7.291906E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(2).real() + 1.224040E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(2).imag() - 4.596569E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).real() - 4.497446E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(2).imag() + 4.072999E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(3).real() + 4.579490E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(3).imag() + 8.590486E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() - 4.313394E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(3).imag() - 4.969719E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(4).real() + 3.356286E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(4).imag() - 3.125121E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(4).real() - 3.171910E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(4).imag() + 3.129068E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(5).real() + 3.149584E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(5).imag() - 3.270358E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(5).real() - 3.105243E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(5).imag() + 3.269355E+03) < 1e-3);

    REQUIRE(
        fabs(result.values.S1(6).real() - 2.25248071E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(6).imag() + 3924.46733361) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).real() + 2.25248071E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).imag() - 3924.46733361) <
        1e-3); // different values using miepython with higher number of terms

    // test single non-magnetic
    refractive_index = std::complex<double>(1.5, -0.5);
    auto size_param_25 = Eigen::VectorXd::LinSpaced(1, 2.5, 2.5);
    result = mie.calculate(size_param_25, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.097071819088392) < 1e-7);
    REQUIRE(fabs(result.values.Qext(0) - 2.562873497454734) < 1e-7);
}

TEST_CASE("LinMie Qext Qsca test small (miepython)", "[sasktran2][mie]") {
    // MieV0 test case 5
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param_99 = Eigen::VectorXd::LinSpaced(1, 0.099, 0.099);
    auto refractive_index = std::complex<double>(0.75, 0);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);

    auto result =
        mie.calculate(size_param_99, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.000007) < 1e-6);

    // MIEV0 test case 6
    auto size_param_101 = Eigen::VectorXd::LinSpaced(1, 0.101, 0.101);
    result = mie.calculate(size_param_101, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.000008) < 1e-6);

    refractive_index = std::complex<double>(1.5, -1);
    auto size_param_55 = Eigen::VectorXd::LinSpaced(1, 0.055, 0.055);
    result = mie.calculate(size_param_55, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.101491) < 1e-6);

    auto size_param_56 = Eigen::VectorXd::LinSpaced(1, 0.056, 0.056);
    result = mie.calculate(size_param_56, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.103347) < 1e-6);

    refractive_index = std::complex<double>(1e-10, -1e-10);
    result = mie.calculate(size_param_99, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 6.329394366790141e-05) < 1e-6);

    result = mie.calculate(size_param_101, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 6.853305788484561e-05) < 1e-6);
}

TEST_CASE("LinMie S1 S2 (miepython)", "[sasktran2][mie]") {
    // MieV0 test case 5
    auto mie = sasktran2::mie::LinearizedMie();
    // auto size_param_1 = Eigen::VectorXd::LinSpaced(1, 1.0, 1.0);
    auto size_param_1 = Eigen::VectorXd::LinSpaced(1, 1.0, 1.0);
    auto refractive_index = std::complex<double>(1.5, -1.0);
    auto angles = Eigen::VectorXd::LinSpaced(7, 0.0, 180.0);
    auto cos_angles = cos(angles.array() / 180.0 * ((double)M_PI));

    auto result =
        mie.calculate(size_param_1, refractive_index, cos_angles, true);
    auto S1 = result.values.S1;
    auto S2 = result.values.S2;
    // auto S1 = result.values.S1.array()*sqrt(((double)
    // M_PI)*size_param_1(0)*size_param_1(0)*result.values.Qext(0)); auto S2 =
    // result.values.S2.array()*sqrt(((double)
    // M_PI)*size_param_1(0)*size_param_1(0)*result.values.Qext(0));

    REQUIRE(fabs(S1(0).real() - 0.584080) < 1e-6);
    REQUIRE(fabs(S1(0).imag() - 0.190515) < 1e-6);
    REQUIRE(fabs(S2(0).real() - 0.584080) < 1e-6);
    REQUIRE(fabs(S2(0).imag() - 0.190515) < 1e-6);

    REQUIRE(fabs(S1(1).real() - 0.565702) < 1e-6);
    REQUIRE(fabs(S1(1).imag() - 0.187200) < 1e-6);
    REQUIRE(fabs(S2(1).real() - 0.500161) < 1e-6);
    REQUIRE(fabs(S2(1).imag() - 0.145611) < 1e-6);

    REQUIRE(fabs(S1(2).real() - 0.517525) < 1e-6);
    REQUIRE(fabs(S1(2).imag() - 0.178443) < 1e-6);
    REQUIRE(fabs(S2(2).real() - 0.287964) < 1e-6);
    REQUIRE(fabs(S2(2).imag() - 0.041054) < 1e-6);

    REQUIRE(fabs(S1(3).real() - 0.456340) < 1e-6);
    REQUIRE(fabs(S1(3).imag() - 0.167167) < 1e-6);
    REQUIRE(fabs(S2(3).real() - 0.0362285) < 1e-6);
    REQUIRE(fabs(S2(3).imag() + 0.0618265) < 1e-6);

    REQUIRE(fabs(S1(4).real() - 0.400212) < 1e-6);
    REQUIRE(fabs(S1(4).imag() - 0.156643) < 1e-6);
    REQUIRE(fabs(S2(4).real() + 0.174875) < 1e-6);
    REQUIRE(fabs(S2(4).imag() + 0.122959) < 1e-6);

    REQUIRE(fabs(S1(5).real() - 0.362157) < 1e-6);
    REQUIRE(fabs(S1(5).imag() - 0.149391) < 1e-6);
    REQUIRE(fabs(S2(5).real() + 0.305682) < 1e-6);
    REQUIRE(fabs(S2(5).imag() + 0.143846) < 1e-6);

    REQUIRE(fabs(S1(6).real() - 0.348844) < 1e-6);
    REQUIRE(fabs(S1(6).imag() - 0.146829) < 1e-6);
    REQUIRE(fabs(S2(6).real() + 0.348844) < 1e-6);
    REQUIRE(fabs(S2(6).imag() + 0.146829) < 1e-6);

    auto size_param_7086 = Eigen::VectorXd::LinSpaced(1, 0.7086, 0.7086);
    refractive_index = std::complex<double>(1.507, -0.002);
    auto cos_angles_2 = Eigen::VectorXd::LinSpaced(1, -1, -1);

    result =
        mie.calculate(size_param_7086, refractive_index, cos_angles_2, true);
    auto S1_2 = result.values.S1.array() /
                sqrt(((double)M_PI) * size_param_7086(0) * size_param_7086(0) *
                     result.values.Qext(0));
    auto S2_2 = result.values.S2.array() /
                sqrt(((double)M_PI) * size_param_7086(0) * size_param_7086(0) *
                     result.values.Qext(0));

    REQUIRE(fabs(S1_2(0).real() - 0.02452300864212876) < 1e-6);
    REQUIRE(fabs(S1_2(0).imag() - 0.29539154027629805) < 1e-6);
    REQUIRE(fabs(S2_2(0).real() + 0.02452300864212876) < 1e-6);
    REQUIRE(fabs(S2_2(0).imag() + 0.29539154027629805) < 1e-6);
}

TEST_CASE("Benchmark for 1000 x cases", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    auto refractive_index = std::complex<double>(1.55, 0);
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);
    // auto size_param = Eigen::VectorXd::LinSpaced(1000, 0.1, 10000);
    auto size_param = Eigen::VectorXd::LinSpaced(1000, 0.1, 10);

    BENCHMARK("Test") {
        return mie.calculate(size_param, refractive_index, cos_angles, true);
    };
}
