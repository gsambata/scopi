#include <xtensor/xmath.hpp>
#include <scopi/object/superellipsoid.hpp>
#include "../mosek_solver.hpp"
#include <random>

int main()
{
    constexpr std::size_t dim = 2;
    double PI = xt::numeric_constants<double>::PI;
    double dt = .05;
    std::size_t total_it = 1000;
    scopi::scopi_container<dim> particles;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distrib_e(1.,1.);
    std::uniform_real_distribution<double> distrib_r(0.1,1.0);
    std::uniform_real_distribution<double> distrib_gp1(2.0,10.0);
    std::uniform_real_distribution<double> distrib_gp2(-10.0,-2.0);
    std::uniform_real_distribution<double> distrib_y(-10.0,10.0);
    std::uniform_real_distribution<double> distrib_rot(0,PI);

    // scopi::superellipsoid<dim> s0({{0.0, 0.}}, {scopi::quaternion(-PI/4)}, {{.01, .01}}, {{0.2}});
    // s0.print();
    // scopi::superellipsoid<dim> s1({{-0.2, 0.}}, {scopi::quaternion(PI/4)}, {{.1, .05}}, {{0.2}});
    // s1.print();
    // scopi::superellipsoid<dim> s2({{0.2, 0.}}, {scopi::quaternion(-PI/4)}, {{.1, .05}}, {{1.}});
    // s2.print();
    // particles.push_back(s0, {{0, 0}}, {{0., 0}}, 0, 0, {{0, 0}});
    // particles.push_back(s1, {{0, 0}}, {{0.25, 0}}, 0, 0, {{0, 0}});
    // particles.push_back(s2, {{0, 0}}, {{-0.25, 0}}, 0, 0, {{0, 0}});

    int n = 10;

    for (int i=0;i<n;++i){

      const double e = distrib_e(generator);

      std::cout << "i = " << i << " e = " << e << std::endl;

      scopi::superellipsoid<dim> s1({{distrib_gp1(generator), distrib_y(generator)}},
        {scopi::quaternion(distrib_rot(generator))}, {{distrib_r(generator), distrib_r(generator)}},
        {{e}});
      // s1.print();
      particles.push_back(s1, {{0, 0}}, {{-0.25, 0}}, 0, 0, {{0, 0}});

      scopi::superellipsoid<dim> s2({{distrib_gp2(generator), distrib_y(generator)}},
        {scopi::quaternion(distrib_rot(generator))}, {{distrib_r(generator), distrib_r(generator)}},
        {{e}});
      // s2.print();
      particles.push_back(s2, {{0, 0}}, {{0.25, 0}}, 0, 0, {{0, 0}});

    }

    // std::cout << "pos = " << particles.pos() << std::endl << std::endl;
    // std::cout << "vd = " << particles.vd() << std::endl << std::endl;
    // std::cout << "q = " << particles.q() << std::endl << std::endl;

    // exit(0);



    // scopi::superellipsoid<dim> s0({{-3.453169298559871,-6.376584982229874}}, {scopi::quaternion(0.9148915951510761)}, {{0.6286366013308933,0.8170145642109674}}, {{0.20927894118442245}});
    // particles.push_back(s0, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s1({{-11.983853712122265,-10.69108824699438}}, {scopi::quaternion(2.5549070456142684)}, {{0.802896838549845,0.4377877124900434}}, {{0.5634513805468084}});
    // particles.push_back(s1, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s2({{-8.044979640195024,-14.666350085848725}}, {scopi::quaternion(0.22716641740589902)}, {{0.6312633837621274,0.8885636937819419}}, {{0.24624446234199635}});
    // particles.push_back(s2, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s3({{-2.52460377004906,-2.4482835671467242}}, {scopi::quaternion(1.2766960641591898)}, {{0.9071125359753474,0.8099593398602669}}, {{0.4113091081381639}});
    // particles.push_back(s3, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s4({{-8.07450573750073,-2.102530116850602}}, {scopi::quaternion(0.7750850646714099)}, {{0.7072913438190973,0.9974345679628966}}, {{0.7966713836646275}});
    // particles.push_back(s4, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s5({{-10.040483216398012,7.421037391527122}}, {scopi::quaternion(0.2823286251978471)}, {{0.9665465149938313,0.3862999101660598}}, {{0.3811291581210021}});
    // particles.push_back(s5, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s6({{-8.05782527022027,1.283177541186685}}, {scopi::quaternion(1.9368314576146641)}, {{0.5324835436027966,0.6138393019914956}}, {{0.6909259148910369}});
    // particles.push_back(s6, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s7({{-13.191212153264171,-4.977340143756619}}, {scopi::quaternion(1.4193422148464592)}, {{0.31313456507023385,0.5807591867264993}}, {{0.27134177703967555}});
    // particles.push_back(s7, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s8({{-13.943628450348953,-11.39932652194295}}, {scopi::quaternion(2.2707299009514763)}, {{0.5069840589411035,0.9717080765702446}}, {{0.3102208024935312}});
    // particles.push_back(s8, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s9({{-12.083352998895096,-15.448184858912137}}, {scopi::quaternion(2.7383900969360426)}, {{0.826472382608898,0.44040857104211795}}, {{0.5537681834319832}});
    // particles.push_back(s9, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s10({{-4.01705613937596,-15.20189658860765}}, {scopi::quaternion(2.725081108462944)}, {{0.632689170037031,0.7057892951784259}}, {{0.7447859221595212}});
    // particles.push_back(s10, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s11({{-2.9957182036735155,-6.80156974759103}}, {scopi::quaternion(3.0217452180485695)}, {{0.9477704671610745,0.30843031221282435}}, {{0.8026036004332282}});
    // particles.push_back(s11, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s12({{-4.0433565875665884,-17.409559881639897}}, {scopi::quaternion(2.2069623407577676)}, {{0.9671165606598495,0.9105700903200216}}, {{0.8026250678013722}});
    // particles.push_back(s12, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s13({{-9.551747567053882,5.801773809194124}}, {scopi::quaternion(0.5136899950631273)}, {{0.7978916288315255,0.36694130953852017}}, {{0.8172323067136009}});
    // particles.push_back(s13, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s14({{-14.438772516660089,-19.04418783732104}}, {scopi::quaternion(0.22898329384959212)}, {{0.3012844277824087,0.7747000006652776}}, {{0.5433308631382429}});
    // particles.push_back(s14, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s15({{-5.327150334801344,-4.009914793091603}}, {scopi::quaternion(0.7974968289761144)}, {{0.6631906177298053,0.9040051061775365}}, {{0.5498858886321767}});
    // particles.push_back(s15, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s16({{-7.655070988556347,7.134938571401523}}, {scopi::quaternion(2.6427848404349468)}, {{0.5261999351906108,0.6644461661588745}}, {{0.8382141846898112}});
    // particles.push_back(s16, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s17({{-3.624949838599541,11.632410518473066}}, {scopi::quaternion(1.2323485759339516)}, {{0.38331031287100514,0.3103760093907067}}, {{0.4444938047901089}});
    // particles.push_back(s17, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s18({{-3.4293864004904115,-9.31463614187229}}, {scopi::quaternion(1.5636654497163167)}, {{0.7373730216541876,0.8184351919733177}}, {{0.6935569621052406}});
    // particles.push_back(s18, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s19({{-13.44367981206335,-18.06309251917312}}, {scopi::quaternion(0.945151725219552)}, {{0.8841466747607329,0.6253510728134344}}, {{0.2950896776064809}});
    // particles.push_back(s19, {{0, 0}}, {{0.5, 0}}, 0, 0, {{0, 0}});
    // // ---------------------------
    // scopi::superellipsoid<dim> s20({{5.029530378075986,-5.7756146712028045}}, {scopi::quaternion(3.0720455159384854)}, {{0.5422962149300444,0.6772151712995982}}, {{0.4562418722443172}});
    // particles.push_back(s20, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s21({{12.228272591038163,4.549591099501448}}, {scopi::quaternion(1.332378996044977)}, {{0.7912041817726029,0.6176092415680858}}, {{0.5241446903362039}});
    // particles.push_back(s21, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s22({{12.748509696249295,-2.0059024253739786}}, {scopi::quaternion(0.03053993413473552)}, {{0.641757555324519,0.8932341050783386}}, {{0.672532458811171}});
    // particles.push_back(s22, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s23({{10.391415923692051,-18.37079574999061}}, {scopi::quaternion(3.080799910249164)}, {{0.7492769871697428,0.45777370568797104}}, {{0.3164107833035701}});
    // particles.push_back(s23, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s24({{13.354158465383925,-3.644382882392687}}, {scopi::quaternion(0.8604880556389541)}, {{0.8533990592499714,0.7661523184975145}}, {{0.9281262882278523}});
    // particles.push_back(s24, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s25({{3.9996365613499196,-17.066620824501637}}, {scopi::quaternion(2.0761511475139236)}, {{0.947229395929275,0.38616393785339787}}, {{0.8061209371134472}});
    // particles.push_back(s25, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s26({{8.933628634299806,-3.929671528228692}}, {scopi::quaternion(1.2253407978519772)}, {{0.9880420178246008,0.9244478049513043}}, {{0.23411938426100526}});
    // particles.push_back(s26, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s27({{6.188110312049356,-7.450375459297245}}, {scopi::quaternion(2.2597158451298247)}, {{0.7475894417477955,0.30728287783771163}}, {{0.4234071954352744}});
    // particles.push_back(s27, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s28({{12.511856192471067,-18.536382771475477}}, {scopi::quaternion(1.6924872067605659)}, {{0.6026646754385917,0.8554347545267862}}, {{0.3765743662980621}});
    // particles.push_back(s28, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s29({{12.428747165577514,-1.0341744367602779}}, {scopi::quaternion(0.7984259142163838)}, {{0.6574170241366936,0.6509728388012103}}, {{0.5915748895716357}});
    // particles.push_back(s29, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s30({{13.141176050042317,-4.945910991427498}}, {scopi::quaternion(1.5263879819038217)}, {{0.7510838743262769,0.9955293095676048}}, {{0.6639049187419357}});
    // particles.push_back(s30, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s31({{14.678477637780649,-3.7037926059404533}}, {scopi::quaternion(0.1711736735658514)}, {{0.5980234242064096,0.5585917675928025}}, {{0.8106352418590015}});
    // particles.push_back(s31, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s32({{11.816378115469183,-14.647951554996364}}, {scopi::quaternion(1.2578815208004759)}, {{0.6357928320038495,0.931444144004909}}, {{0.8558275806814652}});
    // particles.push_back(s32, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s33({{11.22890863598878,-3.375375275802199}}, {scopi::quaternion(0.4461000526779573)}, {{0.7290956687352723,0.6145703692303737}}, {{0.7003178118322819}});
    // particles.push_back(s33, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s34({{3.0215436532170568,-0.7356548743532301}}, {scopi::quaternion(1.5708928457521787)}, {{0.7136682112738463,0.6149607431077839}}, {{0.8321199593921849}});
    // particles.push_back(s34, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s35({{3.6956846808333133,7.678111337885916}}, {scopi::quaternion(0.3507872182039665)}, {{0.7602546475950756,0.9408856918523847}}, {{0.5712654593026503}});
    // particles.push_back(s35, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s36({{3.2815635577518196,-1.9055577815892804}}, {scopi::quaternion(2.81260687553832)}, {{0.5466712597001407,0.8319655401887318}}, {{0.7642360558139427}});
    // particles.push_back(s36, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s37({{13.168523536686678,-17.235630853815852}}, {scopi::quaternion(0.525697248171641)}, {{0.7376517755815697,0.8339585883696679}}, {{0.42437848691131247}});
    // particles.push_back(s37, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s38({{12.568822411359992,4.950748430224397}}, {scopi::quaternion(1.3398592460785392)}, {{0.865635160657374,0.6578213044430823}}, {{0.8768336036288138}});
    // particles.push_back(s38, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});
    // scopi::superellipsoid<dim> s39({{5.6660770013478015,-11.620722509054108}}, {scopi::quaternion(0.7528661770176106)}, {{0.5348415646322534,0.9251702806134072}}, {{0.4818693954212643}});
    // particles.push_back(s39, {{0, 0}}, {{-0.5, 0}}, 0, 0, {{0, 0}});

    // for(std::size_t i = 0; i < particles.size(); ++i)
    // {
    //     std::cout << " i = " << i <<" particles.size() = " << particles.size() << "\n";
    //     particles[i]->print();
    // }
    // exit(0);

    // scopi::superellipsoid<dim> s0({{0.0, 0.}}, {scopi::quaternion(-PI/4)}, {{.01, .01}}, {{0.2}});
    // scopi::superellipsoid<dim> s1({{-0.2, 0.}}, {scopi::quaternion(PI/4)}, {{.1, .05}}, {{0.2}});
    // scopi::superellipsoid<dim> s2({{0.2, 0.}}, {scopi::quaternion(-PI/4)}, {{.1, .05}}, {{1.}});
    // particles.push_back(s0, {{0, 0}}, {{0., 0}}, 0, 0, {{0, 0}});
    // particles.push_back(s1, {{0, 0}}, {{0.25, 0}}, 0, 0, {{0, 0}});
    // particles.push_back(s2, {{0, 0}}, {{-0.25, 0}}, 0, 0, {{0, 0}});

    // std::size_t active_ptr = 1;
    std::size_t active_ptr = 0; // pas d'obstacles

    mosek_solver(particles, dt, total_it, active_ptr);

    return 0;
}
