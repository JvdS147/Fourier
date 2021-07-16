/* *********************************************
Copyright (c) 2013-2021, Cornelis Jan (Jacco) van de Streek
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of my employers nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL CORNELIS JAN VAN DE STREEK BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************* */

#include "Element.h"
#include "BasicMathsFunctions.h"
#include "Utilities.h"

#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream> // For debugging

// Index 0 is D, for all other elements the index is the atomic number

namespace
{

static const double a1[113] =
                   {
                       0.384220004, 0.384220004, 0.768440008, 0.992789984, 2.22744012, 2.03875995, 1.93018997,
                       12.7912998, 2.86039996, 3.30393004, 3.71271992, 5.26399994, 5.59228992, 5.35047007,
                       5.79410982, 6.92073011, 7.06990004, 9.83957005, 16.8752003, 8.11756039, 8.60272026,
                       9.06482029, 9.54969025, 10.0661001, 10.4757004, 11.2518997, 11.9184999, 12.5158005,
                       13.3239002, 13.9351997, 14.6744003, 15.3411999, 15.4378004, 15.4042997, 15.5372000,
                       15.9933996, 16.8493996, 11.4808998, 11.6163998, 19.0566998, 19.2273006, 19.3495998,
                       19.3885002, 19.3596992, 19.4316006, 19.4524002, 19.5123005, 19.5284004, 19.5527992,
                       19.5872002, 19.6527004, 20.0755005, 20.4608002, 20.7492008, 21.6679001, 22.3162994,
                       27.7488995, 33.2108994, 29.4099998, 22.9220009, 23.4069004, 23.8479996, 24.2241993,
                       24.5147991, 24.4004002, 24.3736000, 24.6193008, 24.3162003, 23.8201008, 23.1385994,
                       22.3027992, 21.1865997, 24.6725006, 28.1756992, 31.0935001, 33.2961006, 34.8666992,
                       35.9454002, 36.8101997, 37.3027000, 37.5186005, 37.6946983, 37.7383003, 37.7142982,
                       37.6296997, 37.4971008, 37.3307991, 37.1902008, 36.9819984, 36.8704987, 36.7753983,
                       37.1456985, 37.2807999,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double a2[113] =
                   {
                       0.363469988, 0.363469988, 0.726939976, 0.874019980, 1.55249000, 1.41490996, 1.87811995,
                       3.28546000, 2.52119994, 3.01752996, 3.52630997, 2.17548990, 2.68206000, 2.92451000,
                       3.22390008, 4.14396000, 5.33979988, 7.53180981, 8.32256031, 7.48061991, 7.50768995,
                       7.55526018, 7.60066986, 7.61420012, 7.51401997, 7.36934996, 7.04848003, 6.63641977,
                       6.18745995, 5.84833002, 5.62816000, 5.74149990, 6.00432014, 6.13722992, 5.98288012,
                       6.02439022, 7.19789982, 9.46903992, 9.73009014, 6.50783014, 10.1378002, 10.8737001,
                       11.8308001, 12.8086996, 13.7308998, 14.6844997, 15.3800001, 16.5811005, 17.5716991,
                       18.7168999, 19.5107994, 19.7765999, 20.0335999, 20.5639992, 21.0084991, 21.1791992,
                       21.3777008, 21.7180996, 22.2427998, 22.2518005, 19.7073002, 17.5534992, 15.9132004,
                       14.8057899, 14.0307999, 13.8648996, 14.2735004, 14.9012003, 15.8795996, 17.1707001,
                       18.7201996, 20.1760006, 17.2294998, 14.4287996, 12.5272999, 12.3497000, 11.9524002,
                       11.9980001, 13.0747004, 14.9306002, 17.0352993, 19.7194996, 21.3393993, 22.4542007,
                       23.1322994, 23.5634995, 23.8932991, 24.1306000, 24.2495003, 24.7131004, 25.2506008,
                       25.2998009, 25.6562996,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double a3[113] =
                   {
                       0.138160005, 0.138160005, 0.276309997, 0.842400014, 1.40059996, 1.11609006, 1.57414997,
                       1.76482999, 1.52960002, 1.35754001, 1.19237006, 1.36689997, 1.72235000, 2.27308989,
                       2.42794991, 2.01696992, 2.23629999, 6.07100010, 6.91325998, 1.07795000, 1.75117004,
                       2.05016994, 2.17223001, 2.23551011, 3.50114989, 3.04106998, 3.34326005, 3.57721996,
                       3.74792004, 4.64221001, 3.92540002, 3.10733008, 3.05157995, 3.74678993, 4.83966017,
                       5.51598978, 4.92564011, 9.16981030, 8.68080997, 4.81523991, 2.48177004, 3.47687006,
                       3.75919008, 3.41371989, 4.26536989, 4.50239992, 5.38329983, 4.99149990, 4.47374010,
                       4.02721977, 3.86894989, 4.30389023, 5.38664007, 6.86157990, 8.43381977, 10.7382002,
                       11.0400000, 11.6222000, 11.9818001, 12.2269001, 12.5016003, 12.7323999, 12.9238005,
                       13.0798998, 13.1753998, 13.2510004, 13.3566999, 13.3894997, 13.3937998, 13.3703003,
                       13.3199997, 13.0531998, 12.8069000, 12.6412001, 12.3768997, 11.2819004, 11.1850996,
                       11.2501001, 11.3323002, 10.3424997, 8.51121044, 6.38290024, 5.17527008, 4.84548998,
                       5.59203005, 7.15953016, 9.02221966, 11.5025997, 11.8718996, 12.3888998, 13.0681000,
                       13.7846003, 14.3500996,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double a4[113] =
                   {
                       0.107859999, 0.107859999, 0.215719998, 0.231010005, 0.582899988, 0.732729971, 0.371080011,
                       0.547089994, 0.813499987, 0.836449981, 0.830799997, 1.08859003, 0.730549991, 1.16531003,
                       1.32149005, 1.53859997, 1.51180005, 1.87127995, 2.18514991, 0.972180009, 0.962159991,
                       1.28744996, 1.75437999, 2.23169994, 1.54902005, 2.27702999, 2.27227998, 2.25643992,
                       2.23195004, 1.44753003, 2.16398001, 2.52764010, 2.93571997, 3.01390004, 2.93548989,
                       2.88716006, 2.91605997, 1.42607999, 2.60985994, 2.84786010, 2.42892003, 1.64515996,
                       1.46772003, 1.99925995, 1.28719997, 1.24740005, 0.810150027, 1.21404004, 1.98562002,
                       2.51451993, 3.14763999, 3.44952011, 3.33079004, 2.97588992, 2.62264991, 1.46162999,
                       2.68185997, 3.17238998, 3.19259000, 2.72430992, 2.72849989, 2.72974992, 2.72835994,
                       2.72477007, 3.24471998, 3.24434996, 2.70316005, 2.69308996, 2.68190002, 2.66981006,
                       2.65701008, 3.21190000, 3.55970001, 3.74435997, 3.79137993, 3.72367001, 3.56435990,
                       3.34312010, 2.31420994, 2.01229000, 2.63339996, 3.00959992, 3.71603990, 4.14815998,
                       4.04218006, 3.45923996, 2.77348995, 1.47979999, 2.72428012, 3.26501012, 3.63790989,
                       3.29610991, 3.30732012,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double b1[113] =
                   {
                       10.9070997, 10.9070997, 10.9070997, 4.33978987, 4.964999855E-02, 23.0888004, 12.7187996,
                       2.064000070E-02, 12.7907000, 11.2650995, 3.91090989, 4.02579021, 4.41141987, 3.48664999,
                       2.57103992, 1.83778000, 1.36590004, -5.300000194E-04, -1.456000004E-02, 12.6683998, 10.2636003,
                       8.77431011, 7.60579014, 6.67720985, 6.01658010, 5.34817982, 4.87393999, 4.48994017,
                       4.17742014, 3.97779012, 3.71485996, 3.63867998, 3.39715004, 3.07517004, 2.71530008,
                       2.35650992, 2.01855993, 1.08140004, 1.85573995, 1.24615002, 1.15488005, 1.06625998,
                       0.978770018, 0.893559992, 0.820919991, 0.750190020, 0.685829997, 0.623870015, 0.566039979,
                       0.515100002, 0.466039985, 5.24327993, 4.74224997, 4.27090979, 0.264219999, 0.230920002,
                       0.151519999, 0.110399999, 0.123350002, 2.78604007, 2.71586990, 2.65745997, 2.60993004,
                       2.57255006, 2.47491002, 2.46637011, 2.52207994, 2.52724004, 2.54418993, 2.57319999,
                       2.61392999, 0.886539996, 0.973999977, 1.04033995, 1.07885003, 1.09315002, 1.08840001,
                       1.06923997, 1.04421997, 1.00810003, 0.964550018, 0.922630012, 0.877550006, 0.832220018,
                       0.786400020, 0.740119994, 0.693539977, 0.653029978, 0.603940010, 0.564580023, 0.525099993,
                       0.520200014, 0.502390027,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double b2[113] =
                   {
                       4.30778980, 4.30778980, 4.30778980, 1.26005995, 42.9165001, 0.978479981, 28.6497993,
                       10.7018003, 5.42789984, 4.66504002, 9.63125992, 10.4796000, 1.36548996, 1.20535004,
                       34.1775017, 27.0198002, 19.8278999, 1.11118996, 0.833100021, 0.764090002, 0.627939999,
                       0.533060014, 0.458990008, 0.403219998, 0.374260008, 0.343730003, 0.340229988, 0.354589999,
                       0.386819988, 0.445549995, 0.500329971, 0.656400025, 0.730970025, 0.741129994, 0.689260006,
                       19.7392998, 18.0408993, 18.2800007, 14.6108999, 9.68019009, 10.7876997, 10.5977001,
                       10.0885000, 9.27497005, 8.97737026, 8.42621994, 7.95713997, 7.39504004, 6.79629993,
                       6.29430008, 5.76320982, 0.418579996, 0.370409995, 0.319599986, 3.83525991, 3.49464011,
                       3.09817004, 2.83641005, 2.74836993, 0.180150002, 0.209500000, 0.247799993, 0.294750005,
                       0.349299997, 0.402300000, 0.475169986, 0.545560002, 0.615719974, 0.684449971, 0.749480009,
                       0.808679998, 2.68610001, 2.89037991, 3.20783997, 12.8331003, 13.2559004, 13.8042002,
                       5.43443012, 6.07340002, 6.52549982, 6.65785980, 6.78247976, 6.58964014, 6.27051020,
                       5.86643982, 5.42693996, 4.98695993, 4.61304998, 4.17856979, 3.88775992, 3.61658001,
                       3.66300011, 3.58561993,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double b3[113] =
                   {
                       1.33126998, 1.33126998, 1.33126998, 98.7088013, 1.66378999, 59.8984985, 0.596449971,
                       30.7772999, 0.328200012, 0.337599993, 0.404830009, 0.842220008, 93.4885025, 42.6050987,
                       0.869369984, 0.213180006, 9.210000187E-02, 18.0846004, 14.9176998, 211.222000, 149.300995,
                       123.879997, 109.098999, 98.5953979, 19.0653992, 17.4088993, 15.9329996, 14.8402004,
                       14.0122995, 13.3971004, 12.8862000, 16.0718994, 18.9533005, 21.0014000, 21.0079002,
                       0.581430018, 0.397410005, 2.38825011, 0.898519993, 18.8903008, 120.125999, 32.6174011,
                       31.9738007, 32.3512993, 28.2621002, 26.1564007, 23.1807995, 22.2282009, 21.2907009,
                       22.7308006, 24.0627003, 26.0177994, 27.3458004, 27.3185997, 26.2297001, 25.1863995,
                       20.6774006, 19.3885994, 18.3794003, 17.6662998, 16.9122009, 16.2462997, 15.6554003,
                       15.1280003, 14.4670000, 14.0424004, 13.8486996, 13.5040998, 13.1932001, 12.9125996,
                       12.6590004, 12.2746000, 12.2896996, 12.5053997, 3.63298011, 4.16735983, 4.79179001,
                       14.4982996, 15.7018003, 16.5100002, 16.8437996, 19.2434998, 21.2437000, 24.4692993,
                       27.8677998, 29.8349991, 30.0338001, 29.2597008, 24.3782005, 23.1506004, 22.3409996,
                       20.6539001, 19.6341991,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double b4[113] =
                   {
                       25.6847992, 25.6847992, 25.6847992, 212.087997, 100.361000, 8.538000286E-02, 65.0336990,
                       1.48044002, 34.9406013, 27.9897995, 23.9545994, 133.617004, 32.5280991, 107.169998,
                       85.3410034, 67.1085968, 55.2282982, 45.3666000, 37.2256012, 37.2727013, 60.2274017,
                       36.8889999, 27.5715008, 22.5720005, 97.4598999, 84.2138977, 79.0338974, 74.7351990,
                       71.1194992, 74.1604996, 65.4070969, 70.7609024, 63.7969017, 57.7445984, 52.4308014,
                       47.3322983, 42.5054016, 185.292999, 139.830002, 121.352997, 33.3722000, 120.397003,
                       117.931999, 107.405998, 111.500999, 107.779999, 65.9294968, 100.225998, 85.2777023,
                       88.5674973, 78.1532974, 70.1645966, 65.0572968, 61.5374985, 58.4830017, 232.828995,
                       178.819000, 144.438004, 139.602997, 160.914993, 156.556000, 152.682007, 149.220993,
                       146.102997, 119.738037, 117.445999, 138.384995, 136.246002, 134.281998, 132.468002,
                       130.783005, 107.127998, 93.4381027, 85.0183029, 79.7647018, 76.6561966, 75.1399002,
                       74.7918015, 73.8375015, 76.9116974, 76.7228012, 85.9266968, 78.8094025, 72.1557999,
                       68.1616974, 66.3563995, 65.5799026, 257.964996, 200.024002, 161.725998, 139.164001,
                       150.973007, 146.632996,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double c[113] =
                   {
                       6.250000093E-03, 6.250000093E-03, 1.248999964E-02, 5.987999961E-02, -1.76338995, -0.304089993, 0.246370003,
                       -11.3926001, 0.274199992, 0.483980000, 0.737280011, 1.09912002, 1.26882994, 1.28489006,
                       1.23139000, 0.378699988, -0.159199998, -8.31429958, -16.2971992, 1.35009003, 1.17429996,
                       1.03849006, 0.917620003, 0.845740020, 0.952260017, 1.05194998, 1.40818000, 1.91452003,
                       2.49899006, 3.11685991, 3.59838009, 4.26842022, 4.56067991, 4.69149017, 4.70026016,
                       4.57601976, 4.10864019, 5.43920994, 5.34841013, 5.76120996, 5.71886015, 5.65073013,
                       5.55046988, 5.41555977, 5.28191996, 5.11007023, 4.91426992, 4.68113995, 4.41158009,
                       4.14542007, 3.81226993, 3.38880992, 2.78462005, 1.84739006, 0.266350001, -0.707090020,
                       -6.85854006, -12.7404003, -8.84560013, -1.13929999, 1.64038002, 4.12018013, 6.19355011,
                       7.85730982, 9.12487984, 10.2419996, 11.0290003, 11.6816998, 12.2061996, 12.6322002,
                       12.9818001, 13.3488998, 13.7048998, 13.9823999, 14.1842003, 14.3239002, 14.4097004,
                       14.4448996, 14.4525995, 14.3992004, 14.2910995, 14.1800003, 14.0202999, 13.8301001,
                       13.5991001, 13.3183002, 12.9796000, 12.6868000, 12.1641998, 11.7483997, 11.2496996,
                       11.4561005, 11.3864002,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                        0.0000000
                   };

static const double Van_der_Waals_radii[113] =
                   {
                       1.20, 1.20, 1.40, 1.82, 2.00, 2.00, 1.70, 1.55, 1.52, 1.47, 1.54,
                       2.27, 1.73, 2.00, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.00,
                       2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.40, 1.39,
                       1.87, 2.00, 1.85, 1.90, 1.85, 2.02, 2.00, 2.00, 2.00, 2.00,
                       2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.72, 1.58, 1.93, 2.00,
                       2.00, 2.06, 1.98, 2.16, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
                       2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
                       2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.72, 1.66, 1.55,
                       1.96, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
                       2.00, 1.86, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
                       2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
                       2.00, 2.00
                   };

// What is const here?
static const char* element_symbols_1[113] =
                   {
                       "D",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
                       "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
                       "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                       "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
                       "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
                       "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
                       "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                       "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
                       "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                       "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
                       "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
                       "Rg", "Cn"
                   };

static std::vector< std::string > element_symbols( element_symbols_1, element_symbols_1+113 );

static const double atomic_weights[113] =
                   {
                         2.0,    1.008,   4.003,   6.941,   9.012,  10.811,  12.011,  14.007,  15.999,  18.998,  20.180,
                        22.990, 24.305,  26.982,  28.086,  30.974,  32.064,  35.453,  39.948,  39.098,  40.08,
                        44.956, 47.88,   50.942,  51.996,  54.938,  55.847,  58.933,  58.69,   63.546,  65.38,
                        69.72,  72.61,   74.922,  78.96,   79.904,  83.80,   85.47,   87.62,   88.906,  91.22,
                        92.906, 95.94,   98.0,   101.07,  102.91,  106.4,   107.87,  112.41,  114.82,  118.96,
                       121.75, 127.60,  126.90,  131.29,  132.91,  137.33,  138.91,  140.12,  140.91,  144.24,
                       145.0,  150.36,  151.97,  157.25,  158.93,  162.5,   164.93,  167.26,  168.93,  173.04,
                       174.97, 178.49,  180.95,  183.85,  186.21,  190.2,   192.2,   195.08,  196.97,  200.59,
                       204.38, 207.19,  208.98,  209.0,   210.0,   222.0,   223.0,   226.03,  227.03,  232.04,
                       231.04, 238.03,  237.05,  244.0,   243.0,   247.0,   247.0,   251.0,   252.0,   257.0,
                       258.0,  259.0,   260.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
                         0.0,    0.0
                   };

static const double solid_state_volumes[113] =
                   {
                         5.08,    5.08,   10.0,    22.6,    36.0,    13.24,   13.87,   11.8,    11.39,   11.17,   20.0,
                        26.0,    36.0,    39.6,    37.3,    29.5,    25.2,    25.8,    30.0,    36.0,    45.0,
                        42.0,    27.3,    24.0,    28.1,    31.9,    30.4,    29.4,    26.0,    26.9,    39.0,
                        37.8,    41.6,    36.4,    30.3,    32.7,    40.0,    42.0,    47.0,    44.0,    27.0,
                        37.0,    38.0,    38.0,    37.3,    31.2,    35.0,    35.0,    51.0,    55.0,    52.8,
                        48.0,    46.7,    46.2,    45.0,    46.0,    66.0,    58.0,    54.0,    57.0,    50.0,
                        55.0,    50.0,    53.0,    56.0,    45.0,    50.0,    42.0,    54.0,    49.0,    59.0,
                        35.0,    40.0,    43.0,    38.8,    42.7,    41.9,    34.3,    38.0,    43.0,    38.0,
                        54.0,    52.0,    60.0,    50.0,    55.0,    60.0,    70.0,    60.0,    74.0,    56.0,
                        60.0,    58.0,    45.0,    70.0,    17.0,    70.0,    70.0,    70.0,    70.0,    70.0,
                        70.0,    70.0,    70.0,    70.0,    70.0,    70.0,    70.0,    70.0,    70.0,    70.0,
                        70.0,    70.0
                   };

} // namespace

// ********************************************************************************

Element::Element( const size_t atomic_number ) : id_(atomic_number)
{
    if ( atomic_number > 94 )
        throw std::runtime_error( "Element::Element( const size_t ): atomic number > 94." );
}

// ********************************************************************************

Element::Element( std::string symbol )
{
    // We should probably remove all leading and trailing whitespace
    // Check that length is one or two characters
    if ( ( symbol.length() != 1 ) &&
         ( symbol.length() != 2 ) )
        throw std::runtime_error( "Element::Element( std::string ): length of atomic symbol must be 1 or 2." );
    // Get capitalisation right
    symbol[0] = to_upper( symbol[0] );
    if ( symbol.length() == 2 )
        symbol[1] = to_lower( symbol[1] );
    // Search list
    // Minor problem here if symbol = "H ", which will not match anything
    for ( size_t i(0); i != element_symbols.size(); ++i )
    {
        if ( element_symbols[i] == symbol )
        {
            id_ = i;
            return;
        }
    }
    // If no match, throw an exception
    throw std::runtime_error( "Element::Element( std::string ): could not interpret string " + symbol );
}

// ********************************************************************************

size_t Element::atomic_number() const
{
    return ( id_ == 0 ) ? 1 : id_ ;
}

// ********************************************************************************

std::string Element::symbol() const
{
    return element_symbols[ id_ ];
}

// ********************************************************************************

double Element::atomic_weight() const
{
    return atomic_weights[ id_ ];
}

// ********************************************************************************

double Element::Van_der_Waals_radius() const
{
    return Van_der_Waals_radii[ id_ ];
}

// ********************************************************************************

double Element::solid_state_volume() const
{
    return solid_state_volumes[ id_ ];
}

// ********************************************************************************

double Element::scattering_factor( const double sine_theta_over_lambda ) const
{
    if ( a1[ id_ ] == 0.0 )
        throw std::runtime_error( "Element::scattering_factor(): Scattering factor not yet in table." );
    return a1[ id_ ] * exp(-b1[ id_ ]*square( sine_theta_over_lambda )) +
           a2[ id_ ] * exp(-b2[ id_ ]*square( sine_theta_over_lambda )) +
           a3[ id_ ] * exp(-b3[ id_ ]*square( sine_theta_over_lambda )) +
           a4[ id_ ] * exp(-b4[ id_ ]*square( sine_theta_over_lambda )) +
           c[ id_ ];
}

// ********************************************************************************

//namespace
//{
//
//struct scattering_factor_cache_entry
//{
//    size_t id_;
//    double sine_theta_over_lambda_;
//    double scattering_factor_;
//};
//
//std::vector< scattering_factor_cache_entry > scattering_factor_cache_
//
//} // namespace
//
//double scattering_factor( const Element element, const double sine_theta_over_lambda )
//{
//    if ( a1[ element.id() ] == 0.0000000 )
//        throw std::runtime_error( "Element::scattering_factor(): Scattering factor not yet in table." );
//    // Find in cache
//    // Use of == to test for equality of doubles is deliberate
//    // Manual loop enrolling
//    if ( ( scattering_factor_cache_[0].id_ == id_ ) && ( scattering_factor_cache_[0].sine_theta_over_lambda_ == sine_theta_over_lambda ) )
//        return scattering_factor_cache_[0].scattering_factor_;
//    if ( ( scattering_factor_cache_[1].id_ == id_ ) && ( scattering_factor_cache_[1].sine_theta_over_lambda_ == sine_theta_over_lambda ) )
//        return scattering_factor_cache_[1].scattering_factor_;
//    if ( ( scattering_factor_cache_[2].id_ == id_ ) && ( scattering_factor_cache_[2].sine_theta_over_lambda_ == sine_theta_over_lambda ) )
//        return scattering_factor_cache_[2].scattering_factor_;
//    if ( ( scattering_factor_cache_[3].id_ == id_ ) && ( scattering_factor_cache_[3].sine_theta_over_lambda_ == sine_theta_over_lambda ) )
//        return scattering_factor_cache_[3].scattering_factor_;
//    return a1[ id_ ] * exp(-b1[ id_ ]*square( sine_theta_over_lambda )) +
//           a2[ id_ ] * exp(-b2[ id_ ]*square( sine_theta_over_lambda )) +
//           a3[ id_ ] * exp(-b3[ id_ ]*square( sine_theta_over_lambda )) +
//           a4[ id_ ] * exp(-b4[ id_ ]*square( sine_theta_over_lambda )) +
//           c[ id_ ];
//}

// ********************************************************************************

// The ususal order of elements: C, H, D, rest alphabetical, B before Br
bool elements_less( const Element & lhs, const Element & rhs )
{
    if ( lhs == rhs )
        return false;
    if ( lhs.id() == 6 )
        return true;
    if ( rhs.id() == 6 )
        return false;
    if ( lhs.id() == 1 )
        return true;
    if ( rhs.id() == 1 )
        return false;
    if ( lhs.id() == 0 )
        return true;
    if ( rhs.id() == 0 )
        return false;
    return ( lhs.symbol() < rhs.symbol() );
}

// ********************************************************************************

bool are_bonded( const Element & lhs, const Element & rhs, const double distance2 )
{
    return ( distance2 < square( 0.5 * ( lhs.Van_der_Waals_radius() + rhs.Van_der_Waals_radius() ) ) );
}

// ********************************************************************************

Element element_from_atom_label( std::string label )
{
    label = strip( label );
    if ( label.empty() )
        throw std::runtime_error( "element_from_atom_label(): string is empty." );
    if ( ! std::isalpha( label[0] ) )
        throw std::runtime_error( "element_from_atom_label(): string must start with alphabetic character." );
    if ( ( label.length() == 1 ) || ( ! std::isalpha( label[1] ) ) )
        return Element( label.substr( 0, 1 ) );
    if ( ( to_upper( label[0] ) == 'O' ) && ( to_lower( label[1] ) == 'w' ) )
        return Element( "O" );
    if ( ( label.length() == 2 ) || ( ! std::isalpha( label[2] ) ) )
        return Element( label.substr( 0, 2 ) );
    // We make the best of it, but someone made a booboo
    std::cout << "Warning: element_from_atom_label(): unexpected format: " << label << std::endl;
    return Element( label.substr( 0, 2 ) );
}

// ********************************************************************************

