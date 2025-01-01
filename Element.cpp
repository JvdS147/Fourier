/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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
#include "StringFunctions.h"

#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream> // For debugging

// Index 0 is D, for all other elements the index is the atomic number

namespace
{

static const double a1[113] =
                    {
                         0.38422000400,  0.38422000400,  0.76844000800,  0.99278998400,  2.22744012000,
                         2.03875995000,  1.93018997000, 12.79129980000,  2.86039996000,  3.30393004000,
                         3.71271992000,  5.26399994000,  5.59228992000,  5.35047007000,  5.79410982000,
                         6.92073011000,  7.06990004000,  9.83957005000, 16.87520030000,  8.11756039000,
                         8.60272026000,  9.06482029000,  9.54969025000, 10.06610010000, 10.47570040000,
                        11.25189970000, 11.91849990000, 12.51580050000, 13.32390020000, 13.93519970000,
                        14.67440030000, 15.34119990000, 15.43780040000, 15.40429970000, 15.53720000000,
                        15.99339960000, 16.84939960000, 11.48089980000, 11.61639980000, 19.05669980000,
                        19.22730060000, 19.34959980000, 19.38850020000, 19.35969920000, 19.43160060000,
                        19.45240020000, 19.51230050000, 19.52840040000, 19.55279920000, 19.58720020000,
                        19.65270040000, 20.07550050000, 20.46080020000, 20.74920080000, 21.66790010000,
                        22.31629940000, 27.74889950000, 33.21089940000, 29.40999980000, 22.92200090000,
                        23.40690040000, 23.84799960000, 24.22419930000, 24.51479910000, 24.40040020000,
                        24.37360000000, 24.61930080000, 24.31620030000, 23.82010080000, 23.13859940000,
                        22.30279920000, 21.18659970000, 24.67250060000, 28.17569920000, 31.09350010000,
                        33.29610060000, 34.86669920000, 35.94540020000, 36.81019970000, 37.30270000000,
                        37.51860050000, 37.69469830000, 37.73830030000, 37.71429820000, 37.62969970000,
                        37.49710080000, 37.33079910000, 37.19020080000, 36.98199840000, 36.87049870000,
                        36.77539830000, 37.14569850000, 37.28079990000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000
                    };

static const double a2[113] =
                    {
                         0.36346998800,  0.36346998800,  0.72693997600,  0.87401998000,  1.55249000000,
                         1.41490996000,  1.87811995000,  3.28546000000,  2.52119994000,  3.01752996000,
                         3.52630997000,  2.17548990000,  2.68206000000,  2.92451000000,  3.22390008000,
                         4.14396000000,  5.33979988000,  7.53180981000,  8.32256031000,  7.48061991000,
                         7.50768995000,  7.55526018000,  7.60066986000,  7.61420012000,  7.51401997000,
                         7.36934996000,  7.04848003000,  6.63641977000,  6.18745995000,  5.84833002000,
                         5.62816000000,  5.74149990000,  6.00432014000,  6.13722992000,  5.98288012000,
                         6.02439022000,  7.19789982000,  9.46903992000,  9.73009014000,  6.50783014000,
                        10.13780020000, 10.87370010000, 11.83080010000, 12.80869960000, 13.73089980000,
                        14.68449970000, 15.38000010000, 16.58110050000, 17.57169910000, 18.71689990000,
                        19.51079940000, 19.77659990000, 20.03359990000, 20.56399920000, 21.00849910000,
                        21.17919920000, 21.37770080000, 21.71809960000, 22.24279980000, 22.25180050000,
                        19.70730020000, 17.55349920000, 15.91320040000, 14.80578990000, 14.03079990000,
                        13.86489960000, 14.27350040000, 14.90120030000, 15.87959960000, 17.17070010000,
                        18.72019960000, 20.17600060000, 17.22949980000, 14.42879960000, 12.52729990000,
                        12.34970000000, 11.95240020000, 11.99800010000, 13.07470040000, 14.93060020000,
                        17.03529930000, 19.71949960000, 21.33939930000, 22.45420070000, 23.13229940000,
                        23.56349950000, 23.89329910000, 24.13060000000, 24.24950030000, 24.71310040000,
                        25.25060080000, 25.29980090000, 25.65629960000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000
                    };

static const double a3[113] =
                    {
                         0.13816000500,  0.13816000500,  0.27630999700,  0.84240001400,  1.40059996000,
                         1.11609006000,  1.57414997000,  1.76482999000,  1.52960002000,  1.35754001000,
                         1.19237006000,  1.36689997000,  1.72235000000,  2.27308989000,  2.42794991000,
                         2.01696992000,  2.23629999000,  6.07100010000,  6.91325998000,  1.07795000000,
                         1.75117004000,  2.05016994000,  2.17223001000,  2.23551011000,  3.50114989000,
                         3.04106998000,  3.34326005000,  3.57721996000,  3.74792004000,  4.64221001000,
                         3.92540002000,  3.10733008000,  3.05157995000,  3.74678993000,  4.83966017000,
                         5.51598978000,  4.92564011000,  9.16981030000,  8.68080997000,  4.81523991000,
                         2.48177004000,  3.47687006000,  3.75919008000,  3.41371989000,  4.26536989000,
                         4.50239992000,  5.38329983000,  4.99149990000,  4.47374010000,  4.02721977000,
                         3.86894989000,  4.30389023000,  5.38664007000,  6.86157990000,  8.43381977000,
                        10.73820020000, 11.04000000000, 11.62220000000, 11.98180010000, 12.22690010000,
                        12.50160030000, 12.73239990000, 12.92380050000, 13.07989980000, 13.17539980000,
                        13.25100040000, 13.35669990000, 13.38949970000, 13.39379980000, 13.37030030000,
                        13.31999970000, 13.05319980000, 12.80690000000, 12.64120010000, 12.37689970000,
                        11.28190040000, 11.18509960000, 11.25010010000, 11.33230020000, 10.34249970000,
                         8.51121044000,  6.38290024000,  5.17527008000,  4.84548998000,  5.59203005000,
                         7.15953016000,  9.02221966000, 11.50259970000, 11.87189960000, 12.38889980000,
                        13.06810000000, 13.78460030000, 14.35009960000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000
                    };

static const double a4[113] =
                    {
                        0.10785999900, 0.10785999900, 0.21571999800, 0.23101000500, 0.58289998800,
                        0.73272997100, 0.37108001100, 0.54708999400, 0.81349998700, 0.83644998100,
                        0.83079999700, 1.08859003000, 0.73054999100, 1.16531003000, 1.32149005000,
                        1.53859997000, 1.51180005000, 1.87127995000, 2.18514991000, 0.97218000900,
                        0.96215999100, 1.28744996000, 1.75437999000, 2.23169994000, 1.54902005000,
                        2.27702999000, 2.27227998000, 2.25643992000, 2.23195004000, 1.44753003000,
                        2.16398001000, 2.52764010000, 2.93571997000, 3.01390004000, 2.93548989000,
                        2.88716006000, 2.91605997000, 1.42607999000, 2.60985994000, 2.84786010000,
                        2.42892003000, 1.64515996000, 1.46772003000, 1.99925995000, 1.28719997000,
                        1.24740005000, 0.81015002700, 1.21404004000, 1.98562002000, 2.51451993000,
                        3.14763999000, 3.44952011000, 3.33079004000, 2.97588992000, 2.62264991000,
                        1.46162999000, 2.68185997000, 3.17238998000, 3.19259000000, 2.72430992000,
                        2.72849989000, 2.72974992000, 2.72835994000, 2.72477007000, 3.24471998000,
                        3.24434996000, 2.70316005000, 2.69308996000, 2.68190002000, 2.66981006000,
                        2.65701008000, 3.21190000000, 3.55970001000, 3.74435997000, 3.79137993000,
                        3.72367001000, 3.56435990000, 3.34312010000, 2.31420994000, 2.01229000000,
                        2.63339996000, 3.00959992000, 3.71603990000, 4.14815998000, 4.04218006000,
                        3.45923996000, 2.77348995000, 1.47979999000, 2.72428012000, 3.26501012000,
                        3.63790989000, 3.29610991000, 3.30732012000, 0.00000000000, 0.00000000000,
                        0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000,
                        0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000,
                        0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000,
                        0.00000000000, 0.00000000000, 0.00000000000
                    };

static const double b1[113] =
                    {
                        10.90709970000, 10.90709970000, 10.90709970000,  4.33978987000,  0.04964999855,
                        23.08880040000, 12.71879960000,  0.02064000070, 12.79070000000, 11.26509950000,
                         3.91090989000,  4.02579021000,  4.41141987000,  3.48664999000,  2.57103992000,
                         1.83778000000,  1.36590004000, -0.00053000002, -0.01456000004, 12.66839980000,
                        10.26360030000,  8.77431011000,  7.60579014000,  6.67720985000,  6.01658010000,
                         5.34817982000,  4.87393999000,  4.48994017000,  4.17742014000,  3.97779012000,
                         3.71485996000,  3.63867998000,  3.39715004000,  3.07517004000,  2.71530008000,
                         2.35650992000,  2.01855993000,  1.08140004000,  1.85573995000,  1.24615002000,
                         1.15488005000,  1.06625998000,  0.97877001800,  0.89355999200,  0.82091999100,
                         0.75019002000,  0.68582999700,  0.62387001500,  0.56603997900,  0.51510000200,
                         0.46603998500,  5.24327993000,  4.74224997000,  4.27090979000,  0.26421999900,
                         0.23092000200,  0.15151999900,  0.11039999900,  0.12335000200,  2.78604007000,
                         2.71586990000,  2.65745997000,  2.60993004000,  2.57255006000,  2.47491002000,
                         2.46637011000,  2.52207994000,  2.52724004000,  2.54418993000,  2.57319999000,
                         2.61392999000,  0.88653999600,  0.97399997700,  1.04033995000,  1.07885003000,
                         1.09315002000,  1.08840001000,  1.06923997000,  1.04421997000,  1.00810003000,
                         0.96455001800,  0.92263001200,  0.87755000600,  0.83222001800,  0.78640002000,
                         0.74011999400,  0.69353997700,  0.65302997800,  0.60394001000,  0.56458002300,
                         0.52509999300,  0.52020001400,  0.50239002700,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000
                    };

static const double b2[113] =
                    {
                         4.30778980000,  4.30778980000,  4.30778980000,  1.26005995000, 42.91650010000,
                         0.97847998100, 28.64979930000, 10.70180030000,  5.42789984000,  4.66504002000,
                         9.63125992000, 10.47960000000,  1.36548996000,  1.20535004000, 34.17750170000,
                        27.01980020000, 19.82789990000,  1.11118996000,  0.83310002100,  0.76409000200,
                         0.62793999900,  0.53306001400,  0.45899000800,  0.40321999800,  0.37426000800,
                         0.34373000300,  0.34022998800,  0.35458999900,  0.38681998800,  0.44554999500,
                         0.50032997100,  0.65640002500,  0.73097002500,  0.74112999400,  0.68926000600,
                        19.73929980000, 18.04089930000, 18.28000070000, 14.61089990000,  9.68019009000,
                        10.78769970000, 10.59770010000, 10.08850000000,  9.27497005000,  8.97737026000,
                         8.42621994000,  7.95713997000,  7.39504004000,  6.79629993000,  6.29430008000,
                         5.76320982000,  0.41857999600,  0.37040999500,  0.31959998600,  3.83525991000,
                         3.49464011000,  3.09817004000,  2.83641005000,  2.74836993000,  0.18015000200,
                         0.20950000000,  0.24779999300,  0.29475000500,  0.34929999700,  0.40230000000,
                         0.47516998600,  0.54556000200,  0.61571997400,  0.68444997100,  0.74948000900,
                         0.80867999800,  2.68610001000,  2.89037991000,  3.20783997000, 12.83310030000,
                        13.25590040000, 13.80420020000,  5.43443012000,  6.07340002000,  6.52549982000,
                         6.65785980000,  6.78247976000,  6.58964014000,  6.27051020000,  5.86643982000,
                         5.42693996000,  4.98695993000,  4.61304998000,  4.17856979000,  3.88775992000,
                         3.61658001000,  3.66300011000,  3.58561993000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,
                         0.00000000000,  0.00000000000,  0.00000000000
                    };

static const double b3[113] =
                    {
                          1.33126998000,   1.33126998000,   1.33126998000,  98.70880130000,   1.66378999000,
                         59.89849850000,   0.59644997100,  30.77729990000,   0.32820001200,   0.33759999300,
                          0.40483000900,   0.84222000800,  93.48850250000,  42.60509870000,   0.86936998400,
                          0.21318000600,   0.09210000187,  18.08460040000,  14.91769980000, 211.22200000000,
                        149.30099500000, 123.87999700000, 109.09899900000,  98.59539790000,  19.06539920000,
                         17.40889930000,  15.93299960000,  14.84020040000,  14.01229950000,  13.39710040000,
                         12.88620000000,  16.07189940000,  18.95330050000,  21.00140000000,  21.00790020000,
                          0.58143001800,   0.39741000500,   2.38825011000,   0.89851999300,  18.89030080000,
                        120.12599900000,  32.61740110000,  31.97380070000,  32.35129930000,  28.26210020000,
                         26.15640070000,  23.18079950000,  22.22820090000,  21.29070090000,  22.73080060000,
                         24.06270030000,  26.01779940000,  27.34580040000,  27.31859970000,  26.22970010000,
                         25.18639950000,  20.67740060000,  19.38859940000,  18.37940030000,  17.66629980000,
                         16.91220090000,  16.24629970000,  15.65540030000,  15.12800030000,  14.46700000000,
                         14.04240040000,  13.84869960000,  13.50409980000,  13.19320010000,  12.91259960000,
                         12.65900040000,  12.27460000000,  12.28969960000,  12.50539970000,   3.63298011000,
                          4.16735983000,   4.79179001000,  14.49829960000,  15.70180030000,  16.51000020000,
                         16.84379960000,  19.24349980000,  21.24370000000,  24.46929930000,  27.86779980000,
                         29.83499910000,  30.03380010000,  29.25970080000,  24.37820050000,  23.15060040000,
                         22.34099960000,  20.65390010000,  19.63419910000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000
                    };

static const double b4[113] =
                    {
                         25.68479920000,  25.68479920000,  25.68479920000, 212.08799700000, 100.36100000000,
                          0.08538000286,  65.03369900000,   1.48044002000,  34.94060130000,  27.98979950000,
                         23.95459940000, 133.61700400000,  32.52809910000, 107.16999800000,  85.34100340000,
                         67.10859680000,  55.22829820000,  45.36660000000,  37.22560120000,  37.27270130000,
                         60.22740170000,  36.88899990000,  27.57150080000,  22.57200050000,  97.45989990000,
                         84.21389770000,  79.03389740000,  74.73519900000,  71.11949920000,  74.16049960000,
                         65.40709690000,  70.76090240000,  63.79690170000,  57.74459840000,  52.43080140000,
                         47.33229830000,  42.50540160000, 185.29299900000, 139.83000200000, 121.35299700000,
                         33.37220000000, 120.39700300000, 117.93199900000, 107.40599800000, 111.50099900000,
                        107.77999900000,  65.92949680000, 100.22599800000,  85.27770230000,  88.56749730000,
                         78.15329740000,  70.16459660000,  65.05729680000,  61.53749850000,  58.48300170000,
                        232.82899500000, 178.81900000000, 144.43800400000, 139.60299700000, 160.91499300000,
                        156.55600000000, 152.68200700000, 149.22099300000, 146.10299700000, 119.73803700000,
                        117.44599900000, 138.38499500000, 136.24600200000, 134.28199800000, 132.46800200000,
                        130.78300500000, 107.12799800000,  93.43810270000,  85.01830290000,  79.76470180000,
                         76.65619660000,  75.13990020000,  74.79180150000,  73.83750150000,  76.91169740000,
                         76.72280120000,  85.92669680000,  78.80940250000,  72.15579990000,  68.16169740000,
                         66.35639950000,  65.57990260000, 257.96499600000, 200.02400200000, 161.72599800000,
                        139.16400100000, 150.97300700000, 146.63299600000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000
                    };

static const double c[113] =
                    {
                          0.00625000009,   0.00625000009,   0.01248999964,   0.05987999961,  -1.76338995000,
                         -0.30408999300,   0.24637000300, -11.39260010000,   0.27419999200,   0.48398000000,
                          0.73728001100,   1.09912002000,   1.26882994000,   1.28489006000,   1.23139000000,
                          0.37869998800,  -0.15919999800,  -8.31429958000, -16.29719920000,   1.35009003000,
                          1.17429996000,   1.03849006000,   0.91762000300,   0.84574002000,   0.95226001700,
                          1.05194998000,   1.40818000000,   1.91452003000,   2.49899006000,   3.11685991000,
                          3.59838009000,   4.26842022000,   4.56067991000,   4.69149017000,   4.70026016000,
                          4.57601976000,   4.10864019000,   5.43920994000,   5.34841013000,   5.76120996000,
                          5.71886015000,   5.65073013000,   5.55046988000,   5.41555977000,   5.28191996000,
                          5.11007023000,   4.91426992000,   4.68113995000,   4.41158009000,   4.14542007000,
                          3.81226993000,   3.38880992000,   2.78462005000,   1.84739006000,   0.26635000100,
                         -0.70709002000,  -6.85854006000, -12.74040030000,  -8.84560013000,  -1.13929999000,
                          1.64038002000,   4.12018013000,   6.19355011000,   7.85730982000,   9.12487984000,
                         10.24199960000,  11.02900030000,  11.68169980000,  12.20619960000,  12.63220020000,
                         12.98180010000,  13.34889980000,  13.70489980000,  13.98239990000,  14.18420030000,
                         14.32390020000,  14.40970040000,  14.44489960000,  14.45259950000,  14.39920040000,
                         14.29109950000,  14.18000030000,  14.02029990000,  13.83010010000,  13.59910010000,
                         13.31830020000,  12.97960000000,  12.68680000000,  12.16419980000,  11.74839970000,
                         11.24969960000,  11.45610050000,  11.38640020000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000,
                          0.00000000000,   0.00000000000,   0.00000000000
                    };

static const double Van_der_Waals_radii[113] =
                    {
                        1.20, 1.20, 1.40, 1.82, 2.00, 2.00, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27,
                        1.73, 2.00, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.00, 2.00, 2.00, 2.00,
                        2.00, 2.00, 2.00, 2.00, 1.63, 1.40, 1.39, 1.87, 2.00, 1.85, 1.90, 1.85,
                        2.02, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.72,
                        1.58, 1.93, 2.00, 2.00, 2.06, 1.98, 2.16, 2.00, 2.00, 2.00, 2.00, 2.00,
                        2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
                        2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.72, 1.66, 1.55, 1.96, 2.00, 2.00,
                        2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.86, 2.00, 2.00, 2.00,
                        2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
                        2.00, 2.00, 2.00, 2.00, 2.00
                    };

// What is const here?
static const char* element_symbols_1[113] =
                   {
                       "D" , "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg",
                       "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn",
                       "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
                       "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                       "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                       "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir",
                       "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                       "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                       "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"
                   };

static std::vector< std::string > element_symbols( element_symbols_1, element_symbols_1+113 );

static const double atomic_weights[113] =
                    {
                          2.000,   1.008,   4.003,   6.941,   9.012,  10.811,  12.011,  14.007,
                         15.999,  18.998,  20.180,  22.990,  24.305,  26.982,  28.086,  30.974,
                         32.064,  35.453,  39.948,  39.098,  40.080,  44.956,  47.880,  50.942,
                         51.996,  54.938,  55.847,  58.933,  58.690,  63.546,  65.380,  69.720,
                         72.610,  74.922,  78.960,  79.904,  83.800,  85.470,  87.620,  88.906,
                         91.220,  92.906,  95.940,  98.000, 101.070, 102.910, 106.400, 107.870,
                        112.410, 114.820, 118.960, 121.750, 127.600, 126.900, 131.290, 132.910,
                        137.330, 138.910, 140.120, 140.910, 144.240, 145.000, 150.360, 151.970,
                        157.250, 158.930, 162.500, 164.930, 167.260, 168.930, 173.040, 174.970,
                        178.490, 180.950, 183.850, 186.210, 190.200, 192.200, 195.080, 196.970,
                        200.590, 204.380, 207.190, 208.980, 209.000, 210.000, 222.000, 223.000,
                        226.030, 227.030, 232.040, 231.040, 238.030, 237.050, 244.000, 243.000,
                        247.000, 247.000, 251.000, 252.000, 257.000, 258.000, 259.000, 260.000,
                          0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,
                          0.000
                    };

static const double solid_state_volumes[113] =
                    {
                         5.08,  5.08, 10.00, 22.60, 36.00, 13.24, 13.87, 11.80, 11.39, 11.17, 20.00,
                        26.00, 36.00, 39.60, 37.30, 29.50, 25.20, 25.80, 30.00, 36.00, 45.00, 42.00,
                        27.30, 24.00, 28.10, 31.90, 30.40, 29.40, 26.00, 26.90, 39.00, 37.80, 41.60,
                        36.40, 30.30, 32.70, 40.00, 42.00, 47.00, 44.00, 27.00, 37.00, 38.00, 38.00,
                        37.30, 31.20, 35.00, 35.00, 51.00, 55.00, 52.80, 48.00, 46.70, 46.20, 45.00,
                        46.00, 66.00, 58.00, 54.00, 57.00, 50.00, 55.00, 50.00, 53.00, 56.00, 45.00,
                        50.00, 42.00, 54.00, 49.00, 59.00, 35.00, 40.00, 43.00, 38.80, 42.70, 41.90,
                        34.30, 38.00, 43.00, 38.00, 54.00, 52.00, 60.00, 50.00, 55.00, 60.00, 70.00,
                        60.00, 74.00, 56.00, 60.00, 58.00, 45.00, 70.00, 17.00, 70.00, 70.00, 70.00,
                        70.00, 70.00, 70.00, 70.00, 70.00, 70.00, 70.00, 70.00, 70.00, 70.00, 70.00,
                        70.00, 70.00, 70.00
                    };

} // namespace

// ********************************************************************************

Element::Element( const size_t atomic_number ) : id_(atomic_number)
{
    if ( atomic_number > 112 )
        throw std::runtime_error( "Element::Element( const size_t ): atomic number > 112." );
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
        throw std::runtime_error( "element_from_atom_label(): string >" + label + "< must start with alphabetic character." );
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

