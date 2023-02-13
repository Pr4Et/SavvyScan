// ***********************************************************************
//
// regs.h                                          (c) Spectrum GmbH, 2006
//
// ***********************************************************************
//
// software register and constants definition for all Spectrum drivers. 
// Please stick to the card manual to see which of the inhere defined 
// registers are used on your hardware.
//
// ***********************************************************************



// ***********************************************************************
// macros for kilo, Mega or Giga as standard version or binary (_B) (2^x)
// ***********************************************************************

#define KILO(k)     ((uint64) 1000 * (k))
#define MEGA(m)     ((uint64) 1000 * 1000 * (m))
#define GIGA(g)     ((uint64) 1000 * 1000 * 1000 * (g))
#define KILO_B(k)   ((uint64) 1024 * (k))
#define MEGA_B(m)   ((uint64) 1024 * 1024 * (m))
#define GIGA_B(g)   ((uint64) 1024 * 1024 * 1024 * (g))




// ***********************************************************************
// card types
// ***********************************************************************

#define TYP_PCIDEVICEID             0x00000000l

// ***** Board Types ***************
#define TYP_EVAL                    0x00000010l
#define TYP_RSDLGA                  0x00000014l
#define TYP_GMG                     0x00000018l
#define TYP_VAN8                    0x00000020l
#define TYP_VAC                     0x00000028l

#define TYP_PCIAUTOINSTALL          0x000000FFl

#define TYP_DAP116                  0x00000100l
#define TYP_PAD82                   0x00000200l
#define TYP_PAD82a                  0x00000210l
#define TYP_PAD82b                  0x00000220l
#define TYP_PCI212                  0x00000300l
#define TYP_PAD1232a                0x00000400l
#define TYP_PAD1232b                0x00000410l
#define TYP_PAD1232c                0x00000420l
#define TYP_PAD1616a                0x00000500l
#define TYP_PAD1616b                0x00000510l
#define TYP_PAD1616c                0x00000520l
#define TYP_PAD1616d                0x00000530l
#define TYP_PAD52                   0x00000600l
#define TYP_PAD242                  0x00000700l
#define TYP_PCK400                  0x00000800l
#define TYP_PAD164_2M               0x00000900l
#define TYP_PAD164_5M               0x00000910l
#define TYP_PCI208                  0x00001000l
#define TYP_CPCI208                 0x00001001l
#define TYP_PCI412                  0x00001100l
#define TYP_PCIDIO32                0x00001200l
#define TYP_PCI248                  0x00001300l
#define TYP_PADCO                   0x00001400l
#define TYP_TRS582                  0x00001500l
#define TYP_PCI258                  0x00001600l


// ------ series and familiy identifiers -----
#define TYP_SERIESMASK              0x00FF0000l     // the series (= type of base card), e.g. MI.xxxx
#define TYP_VERSIONMASK             0x0000FFFFl     // the version, e.g. XX.3012
#define TYP_FAMILYMASK              0x0000FF00l     // the family, e.g. XX.30xx
#define TYP_TYPEMASK                0x000000FFl     // the type, e.g. XX.xx12
#define TYP_SPEEDMASK               0x000000F0l     // the speed grade, e.g. XX.xx1x
#define TYP_CHMASK                  0x0000000Fl     // the channel/modules, e.g. XX.xxx2

#define TYP_MISERIES                0x00000000l
#define TYP_MCSERIES                0x00010000l
#define TYP_MXSERIES                0x00020000l
#define TYP_M2ISERIES               0x00030000l
#define TYP_M2IEXPSERIES            0x00040000l
#define TYP_M3ISERIES               0x00050000l
#define TYP_M3IEXPSERIES            0x00060000l
#define TYP_M4IEXPSERIES            0x00070000l
#define TYP_M4XEXPSERIES            0x00080000l
#define TYP_M2PEXPSERIES            0x00090000l



// ----- MI.20xx, MC.20xx, MX.20xx -----
#define TYP_MI2020                  0x00002020l
#define TYP_MI2021                  0x00002021l
#define TYP_MI2025                  0x00002025l
#define TYP_MI2030                  0x00002030l
#define TYP_MI2031                  0x00002031l

#define TYP_M2I2020                 0x00032020l
#define TYP_M2I2021                 0x00032021l
#define TYP_M2I2025                 0x00032025l
#define TYP_M2I2030                 0x00032030l
#define TYP_M2I2031                 0x00032031l

#define TYP_M2I2020EXP              0x00042020l
#define TYP_M2I2021EXP              0x00042021l
#define TYP_M2I2025EXP              0x00042025l
#define TYP_M2I2030EXP              0x00042030l
#define TYP_M2I2031EXP              0x00042031l

#define TYP_MC2020                  0x00012020l
#define TYP_MC2021                  0x00012021l
#define TYP_MC2025                  0x00012025l
#define TYP_MC2030                  0x00012030l
#define TYP_MC2031                  0x00012031l

#define TYP_MX2020                  0x00022020l
#define TYP_MX2025                  0x00022025l
#define TYP_MX2030                  0x00022030l

// ----- M3i.21xx, M3i.21xx-Exp (8 bit A/D) -----
#define TYP_M3I2120                 0x00052120l     // 1x500M
#define TYP_M3I2122                 0x00052122l     // 1x500M & 2x250M
#define TYP_M3I2130                 0x00052130l     // 1x1G
#define TYP_M3I2132                 0x00052132l     // 1x1G & 2x500M

#define TYP_M3I2120EXP              0x00062120l     // 1x500M
#define TYP_M3I2122EXP              0x00062122l     // 1x500M & 2x250M
#define TYP_M3I2130EXP              0x00062130l     // 1x1G
#define TYP_M3I2132EXP              0x00062132l     // 1x1G & 2x500M

// ----- M4i.22xx-x8 (8 bit A/D) -----
#define TYP_M4I22XX_X8              0x00072200l
#define TYP_M4I2210_X8              0x00072210l     // 1x1.25G
#define TYP_M4I2211_X8              0x00072211l     // 2x1.25G
#define TYP_M4I2212_X8              0x00072212l     // 4x1.25G
#define TYP_M4I2220_X8              0x00072220l     // 1x2.5G
#define TYP_M4I2221_X8              0x00072221l     // 2x2.5G
#define TYP_M4I2223_X8              0x00072223l     // 1x2.5G & 2x1.25G
#define TYP_M4I2230_X8              0x00072230l     // 1x5G
#define TYP_M4I2233_X8              0x00072233l     // 1x5G & 2x2.5G 
#define TYP_M4I2234_X8              0x00072234l     // 1x5G & 2x2.5G & 4x1.25G
#define TYP_M4I2280_X8              0x00072280l     // customer specific variant
#define TYP_M4I2281_X8              0x00072281l     // customer specific variant
#define TYP_M4I2283_X8              0x00072283l     // customer specific variant
#define TYP_M4I2290_X8              0x00072290l     // customer specific variant
#define TYP_M4I2293_X8              0x00072293l     // customer specific variant
#define TYP_M4I2294_X8              0x00072294l     // customer specific variant

// ----- M4x.22xx-x8 (8 bit A/D) -----
#define TYP_M4X22XX_X4              0x00082200l
#define TYP_M4X2210_X4              0x00082210l     // 1x1.25G
#define TYP_M4X2211_X4              0x00082211l     // 2x1.25G
#define TYP_M4X2212_X4              0x00082212l     // 4x1.25G
#define TYP_M4X2220_X4              0x00082220l     // 1x2.5G
#define TYP_M4X2221_X4              0x00082221l     // 2x2.5G
#define TYP_M4X2223_X4              0x00082223l     // 1x2.5G & 2x1.25G
#define TYP_M4X2230_X4              0x00082230l     // 1x5G
#define TYP_M4X2233_X4              0x00082233l     // 1x5G & 2x2.5G 
#define TYP_M4X2234_X4              0x00082234l     // 1x5G & 2x2.5G & 4x1.25G

// ----- M4i.23xx-x8 (7 bit A/D) -----
#define TYP_M4I23XX_X8              0x00072300l
#define TYP_M4I2320_X8              0x00072320l     // 1x2.5G
#define TYP_M4I2321_X8              0x00072321l     // 2x2.5G
#define TYP_M4I2323_X8              0x00072323l     // 1x2.5G & 2x1.25G
#define TYP_M4I2330_X8              0x00072330l     // 1x5G
#define TYP_M4I2333_X8              0x00072333l     // 1x5G & 2x2.5G 
#define TYP_M4I2334_X8              0x00072334l     // 1x5G & 2x2.5G & 4x1.25G

// ----- MI.30xx, MC.30xx, MX.30xx -----
#define TYP_MI3010                  0x00003010l
#define TYP_MI3011                  0x00003011l
#define TYP_MI3012                  0x00003012l
#define TYP_MI3013                  0x00003013l
#define TYP_MI3014                  0x00003014l
#define TYP_MI3015                  0x00003015l
#define TYP_MI3016                  0x00003016l
#define TYP_MI3020                  0x00003020l
#define TYP_MI3021                  0x00003021l
#define TYP_MI3022                  0x00003022l
#define TYP_MI3023                  0x00003023l
#define TYP_MI3024                  0x00003024l
#define TYP_MI3025                  0x00003025l
#define TYP_MI3026                  0x00003026l
#define TYP_MI3027                  0x00003027l
#define TYP_MI3031                  0x00003031l
#define TYP_MI3033                  0x00003033l

#define TYP_M2I3010                 0x00033010l
#define TYP_M2I3011                 0x00033011l
#define TYP_M2I3012                 0x00033012l
#define TYP_M2I3013                 0x00033013l
#define TYP_M2I3014                 0x00033014l
#define TYP_M2I3015                 0x00033015l
#define TYP_M2I3016                 0x00033016l
#define TYP_M2I3020                 0x00033020l
#define TYP_M2I3021                 0x00033021l
#define TYP_M2I3022                 0x00033022l
#define TYP_M2I3023                 0x00033023l
#define TYP_M2I3024                 0x00033024l
#define TYP_M2I3025                 0x00033025l
#define TYP_M2I3026                 0x00033026l
#define TYP_M2I3027                 0x00033027l
#define TYP_M2I3031                 0x00033031l
#define TYP_M2I3033                 0x00033033l

#define TYP_M2I3010EXP              0x00043010l
#define TYP_M2I3011EXP              0x00043011l
#define TYP_M2I3012EXP              0x00043012l
#define TYP_M2I3013EXP              0x00043013l
#define TYP_M2I3014EXP              0x00043014l
#define TYP_M2I3015EXP              0x00043015l
#define TYP_M2I3016EXP              0x00043016l
#define TYP_M2I3020EXP              0x00043020l
#define TYP_M2I3021EXP              0x00043021l
#define TYP_M2I3022EXP              0x00043022l
#define TYP_M2I3023EXP              0x00043023l
#define TYP_M2I3024EXP              0x00043024l
#define TYP_M2I3025EXP              0x00043025l
#define TYP_M2I3026EXP              0x00043026l
#define TYP_M2I3027EXP              0x00043027l
#define TYP_M2I3031EXP              0x00043031l
#define TYP_M2I3033EXP              0x00043033l

#define TYP_MC3010                  0x00013010l
#define TYP_MC3011                  0x00013011l
#define TYP_MC3012                  0x00013012l
#define TYP_MC3013                  0x00013013l
#define TYP_MC3014                  0x00013014l
#define TYP_MC3015                  0x00013015l
#define TYP_MC3016                  0x00013016l
#define TYP_MC3020                  0x00013020l
#define TYP_MC3021                  0x00013021l
#define TYP_MC3022                  0x00013022l
#define TYP_MC3023                  0x00013023l
#define TYP_MC3024                  0x00013024l
#define TYP_MC3025                  0x00013025l
#define TYP_MC3026                  0x00013026l
#define TYP_MC3027                  0x00013027l
#define TYP_MC3031                  0x00013031l
#define TYP_MC3033                  0x00013033l

#define TYP_MX3010                  0x00023010l
#define TYP_MX3011                  0x00023011l
#define TYP_MX3012                  0x00023012l
#define TYP_MX3020                  0x00023020l
#define TYP_MX3021                  0x00023021l
#define TYP_MX3022                  0x00023022l
#define TYP_MX3031                  0x00023031l



// ----- MI.31xx, MC.31xx, MX.31xx -----
#define TYP_MI3110                  0x00003110l
#define TYP_MI3111                  0x00003111l
#define TYP_MI3112                  0x00003112l
#define TYP_MI3120                  0x00003120l
#define TYP_MI3121                  0x00003121l
#define TYP_MI3122                  0x00003122l
#define TYP_MI3130                  0x00003130l
#define TYP_MI3131                  0x00003131l
#define TYP_MI3132                  0x00003132l
#define TYP_MI3140                  0x00003140l

#define TYP_M2I3110                 0x00033110l
#define TYP_M2I3111                 0x00033111l
#define TYP_M2I3112                 0x00033112l
#define TYP_M2I3120                 0x00033120l
#define TYP_M2I3121                 0x00033121l
#define TYP_M2I3122                 0x00033122l
#define TYP_M2I3130                 0x00033130l
#define TYP_M2I3131                 0x00033131l
#define TYP_M2I3132                 0x00033132l

#define TYP_M2I3110EXP              0x00043110l
#define TYP_M2I3111EXP              0x00043111l
#define TYP_M2I3112EXP              0x00043112l
#define TYP_M2I3120EXP              0x00043120l
#define TYP_M2I3121EXP              0x00043121l
#define TYP_M2I3122EXP              0x00043122l
#define TYP_M2I3130EXP              0x00043130l
#define TYP_M2I3131EXP              0x00043131l
#define TYP_M2I3132EXP              0x00043132l

#define TYP_MC3110                  0x00013110l
#define TYP_MC3111                  0x00013111l
#define TYP_MC3112                  0x00013112l
#define TYP_MC3120                  0x00013120l
#define TYP_MC3121                  0x00013121l
#define TYP_MC3122                  0x00013122l
#define TYP_MC3130                  0x00013130l
#define TYP_MC3131                  0x00013131l
#define TYP_MC3132                  0x00013132l

#define TYP_MX3110                  0x00023110l
#define TYP_MX3111                  0x00023111l
#define TYP_MX3120                  0x00023120l
#define TYP_MX3121                  0x00023121l
#define TYP_MX3130                  0x00023130l
#define TYP_MX3131                  0x00023131l



// ----- M3i.32xx, M3i.32xx-Exp (12 bit A/D) -----
#define TYP_M3I3220                 0x00053220l     // 1x250M
#define TYP_M3I3221                 0x00053221l     // 2x250M
#define TYP_M3I3240                 0x00053240l     // 1x500M
#define TYP_M3I3242                 0x00053242l     // 1x500M & 2x250M

#define TYP_M3I3220EXP              0x00063220l     // 1x250M
#define TYP_M3I3221EXP              0x00063221l     // 2x250M
#define TYP_M3I3240EXP              0x00063240l     // 1x500M
#define TYP_M3I3242EXP              0x00063242l     // 1x500M & 2x250M



// ----- MI.40xx, MC.40xx, MX.40xx -----
#define TYP_MI4020                  0x00004020l
#define TYP_MI4021                  0x00004021l
#define TYP_MI4022                  0x00004022l
#define TYP_MI4030                  0x00004030l
#define TYP_MI4031                  0x00004031l
#define TYP_MI4032                  0x00004032l

#define TYP_M2I4020                 0x00034020l
#define TYP_M2I4021                 0x00034021l
#define TYP_M2I4022                 0x00034022l
#define TYP_M2I4028                 0x00034028l
#define TYP_M2I4030                 0x00034030l
#define TYP_M2I4031                 0x00034031l
#define TYP_M2I4032                 0x00034032l
#define TYP_M2I4038                 0x00034038l

#define TYP_M2I4020EXP              0x00044020l
#define TYP_M2I4021EXP              0x00044021l
#define TYP_M2I4022EXP              0x00044022l
#define TYP_M2I4028EXP              0x00044028l
#define TYP_M2I4030EXP              0x00044030l
#define TYP_M2I4031EXP              0x00044031l
#define TYP_M2I4032EXP              0x00044032l
#define TYP_M2I4038EXP              0x00044038l

#define TYP_MC4020                  0x00014020l
#define TYP_MC4021                  0x00014021l
#define TYP_MC4022                  0x00014022l
#define TYP_MC4030                  0x00014030l
#define TYP_MC4031                  0x00014031l
#define TYP_MC4032                  0x00014032l

#define TYP_MX4020                  0x00024020l
#define TYP_MX4021                  0x00024021l
#define TYP_MX4030                  0x00024030l
#define TYP_MX4031                  0x00024031l



// ----- M3i.41xx, M3i.41xx-Exp (14 bit A/D) -----
#define TYP_M3I4110                 0x00054110l     // 1x100M
#define TYP_M3I4111                 0x00054111l     // 2x100M
#define TYP_M3I4120                 0x00054120l     // 1x250M
#define TYP_M3I4121                 0x00054121l     // 2x250M
#define TYP_M3I4140                 0x00054140l     // 1x400M
#define TYP_M3I4142                 0x00054142l     // 1x400M & 2x250M

#define TYP_M3I4110EXP              0x00064110l     // 1x100M
#define TYP_M3I4111EXP              0x00064111l     // 2x100M
#define TYP_M3I4120EXP              0x00064120l     // 1x250M
#define TYP_M3I4121EXP              0x00064121l     // 2x250M
#define TYP_M3I4140EXP              0x00064140l     // 1x400M
#define TYP_M3I4142EXP              0x00064142l     // 1x400M & 2x250M

// ----- M4i.44xx-x8 (generic) -----
#define TYP_M4I44XX_X8              0x00074400      // 

#define TYP_M4I4410_X8              0x00074410      // 2x130M 16bit
#define TYP_M4I4411_X8              0x00074411      // 4x130M 16bit
#define TYP_M4I4420_X8              0x00074420      // 2x250M 16bit
#define TYP_M4I4421_X8              0x00074421      // 4x250M 16bit
#define TYP_M4I4450_X8              0x00074450      // 2x500M 14bit
#define TYP_M4I4451_X8              0x00074451      // 4x500M 14bit
#define TYP_M4I4470_X8              0x00074470      // 2x180M 16bit
#define TYP_M4I4471_X8              0x00074471      // 4x180M 16bit
#define TYP_M4I4480_X8              0x00074480      // 2x400M 14bit
#define TYP_M4I4481_X8              0x00074481      // 4x400M 14bit

// ----- M4x.44xx-x4 (14/16 bit A/D) -----
#define TYP_M4X44XX_X4              0x00084400      // generic
#define TYP_M4X4410_X4              0x00084410      // 2x130M 16bit
#define TYP_M4X4411_X4              0x00084411      // 4x130M 16bit
#define TYP_M4X4420_X4              0x00084420      // 2x250M 16bit
#define TYP_M4X4421_X4              0x00084421      // 4x250M 16bit
#define TYP_M4X4450_X4              0x00084450      // 2x500M 14bit
#define TYP_M4X4451_X4              0x00084451      // 4x500M 14bit
#define TYP_M4X4470_X4              0x00084470      // 2x180M 16bit
#define TYP_M4X4471_X4              0x00084471      // 4x180M 16bit
#define TYP_M4X4480_X4              0x00084480      // 2x400M 14bit
#define TYP_M4X4481_X4              0x00084481      // 4x400M 14bit


// ----- MI.45xx, MC.45xx, MX.45xx -----
#define TYP_MI4520                  0x00004520l
#define TYP_MI4521                  0x00004521l
#define TYP_MI4530                  0x00004530l
#define TYP_MI4531                  0x00004531l
#define TYP_MI4540                  0x00004540l
#define TYP_MI4541                  0x00004541l

#define TYP_M2I4520                 0x00034520l
#define TYP_M2I4521                 0x00034521l
#define TYP_M2I4530                 0x00034530l
#define TYP_M2I4531                 0x00034531l
#define TYP_M2I4540                 0x00034540l
#define TYP_M2I4541                 0x00034541l

#define TYP_MC4520                  0x00014520l
#define TYP_MC4521                  0x00014521l
#define TYP_MC4530                  0x00014530l
#define TYP_MC4531                  0x00014531l
#define TYP_MC4540                  0x00014540l
#define TYP_MC4541                  0x00014541l

#define TYP_MX4520                  0x00024520l
#define TYP_MX4530                  0x00024530l
#define TYP_MX4540                  0x00024540l



// ----- MI.46xx, MC.46xx, MX.46xx -----
#define TYP_MI4620                  0x00004620l
#define TYP_MI4621                  0x00004621l
#define TYP_MI4622                  0x00004622l
#define TYP_MI4630                  0x00004630l
#define TYP_MI4631                  0x00004631l
#define TYP_MI4632                  0x00004632l
#define TYP_MI4640                  0x00004640l
#define TYP_MI4641                  0x00004641l
#define TYP_MI4642                  0x00004642l
#define TYP_MI4650                  0x00004650l
#define TYP_MI4651                  0x00004651l
#define TYP_MI4652                  0x00004652l

#define TYP_M2I4620                 0x00034620l
#define TYP_M2I4621                 0x00034621l
#define TYP_M2I4622                 0x00034622l
#define TYP_M2I4630                 0x00034630l
#define TYP_M2I4631                 0x00034631l
#define TYP_M2I4632                 0x00034632l
#define TYP_M2I4640                 0x00034640l
#define TYP_M2I4641                 0x00034641l
#define TYP_M2I4642                 0x00034642l
#define TYP_M2I4650                 0x00034650l
#define TYP_M2I4651                 0x00034651l
#define TYP_M2I4652                 0x00034652l

#define TYP_M2I4620EXP              0x00044620l
#define TYP_M2I4621EXP              0x00044621l
#define TYP_M2I4622EXP              0x00044622l
#define TYP_M2I4630EXP              0x00044630l
#define TYP_M2I4631EXP              0x00044631l
#define TYP_M2I4632EXP              0x00044632l
#define TYP_M2I4640EXP              0x00044640l
#define TYP_M2I4641EXP              0x00044641l
#define TYP_M2I4642EXP              0x00044642l
#define TYP_M2I4650EXP              0x00044650l
#define TYP_M2I4651EXP              0x00044651l
#define TYP_M2I4652EXP              0x00044652l

#define TYP_MC4620                  0x00014620l
#define TYP_MC4621                  0x00014621l
#define TYP_MC4622                  0x00014622l
#define TYP_MC4630                  0x00014630l
#define TYP_MC4631                  0x00014631l
#define TYP_MC4632                  0x00014632l
#define TYP_MC4640                  0x00014640l
#define TYP_MC4641                  0x00014641l
#define TYP_MC4642                  0x00014642l
#define TYP_MC4650                  0x00014650l
#define TYP_MC4651                  0x00014651l
#define TYP_MC4652                  0x00014652l

#define TYP_MX4620                  0x00024620l
#define TYP_MX4621                  0x00024621l
#define TYP_MX4630                  0x00024630l
#define TYP_MX4631                  0x00024631l
#define TYP_MX4640                  0x00024640l
#define TYP_MX4641                  0x00024641l
#define TYP_MX4650                  0x00024650l
#define TYP_MX4651                  0x00024651l



// ----- MI.47xx, MC.47xx, MX.47xx -----
#define TYP_MI4710                  0x00004710l
#define TYP_MI4711                  0x00004711l
#define TYP_MI4720                  0x00004720l
#define TYP_MI4721                  0x00004721l
#define TYP_MI4730                  0x00004730l
#define TYP_MI4731                  0x00004731l
#define TYP_MI4740                  0x00004740l
#define TYP_MI4741                  0x00004741l

#define TYP_M2I4710                 0x00034710l
#define TYP_M2I4711                 0x00034711l
#define TYP_M2I4720                 0x00034720l
#define TYP_M2I4721                 0x00034721l
#define TYP_M2I4730                 0x00034730l
#define TYP_M2I4731                 0x00034731l
#define TYP_M2I4740                 0x00034740l
#define TYP_M2I4741                 0x00034741l

#define TYP_M2I4710EXP              0x00044710l
#define TYP_M2I4711EXP              0x00044711l
#define TYP_M2I4720EXP              0x00044720l
#define TYP_M2I4721EXP              0x00044721l
#define TYP_M2I4730EXP              0x00044730l
#define TYP_M2I4731EXP              0x00044731l
#define TYP_M2I4740EXP              0x00044740l
#define TYP_M2I4741EXP              0x00044741l

#define TYP_MC4710                  0x00014710l
#define TYP_MC4711                  0x00014711l
#define TYP_MC4720                  0x00014720l
#define TYP_MC4721                  0x00014721l
#define TYP_MC4730                  0x00014730l
#define TYP_MC4731                  0x00014731l

#define TYP_MX4710                  0x00024710l
#define TYP_MX4720                  0x00024720l
#define TYP_MX4730                  0x00024730l



// ----- M3i.48xx, M3i.48xx-Exp (16 bit A/D) -----
#define TYP_M3I4830                 0x00054830l     
#define TYP_M3I4831                 0x00054831l    
#define TYP_M3I4840                 0x00054840l     
#define TYP_M3I4841                 0x00054841l    
#define TYP_M3I4860                 0x00054860l     
#define TYP_M3I4861                 0x00054861l    

#define TYP_M3I4830EXP              0x00064830l     
#define TYP_M3I4831EXP              0x00064831l    
#define TYP_M3I4840EXP              0x00064840l     
#define TYP_M3I4841EXP              0x00064841l    
#define TYP_M3I4860EXP              0x00064860l     
#define TYP_M3I4861EXP              0x00064861l    



// ----- MI.46xx, MC.46xx, MX.46xx -----
#define TYP_MI4911                  0x00004911l
#define TYP_MI4912                  0x00004912l
#define TYP_MI4931                  0x00004931l
#define TYP_MI4932                  0x00004932l
#define TYP_MI4960                  0x00004960l
#define TYP_MI4961                  0x00004961l
#define TYP_MI4963                  0x00004963l
#define TYP_MI4964                  0x00004964l

#define TYP_MC4911                  0x00014911l
#define TYP_MC4912                  0x00014912l
#define TYP_MC4931                  0x00014931l
#define TYP_MC4932                  0x00014932l
#define TYP_MC4960                  0x00014960l
#define TYP_MC4961                  0x00014961l
#define TYP_MC4963                  0x00014963l
#define TYP_MC4964                  0x00014964l

#define TYP_MX4911                  0x00024911l
#define TYP_MX4931                  0x00024931l
#define TYP_MX4960                  0x00024960l
#define TYP_MX4963                  0x00024963l

#define TYP_M2I4911                 0x00034911l
#define TYP_M2I4912                 0x00034912l
#define TYP_M2I4931                 0x00034931l
#define TYP_M2I4932                 0x00034932l
#define TYP_M2I4960                 0x00034960l
#define TYP_M2I4961                 0x00034961l
#define TYP_M2I4963                 0x00034963l
#define TYP_M2I4964                 0x00034964l

#define TYP_M2I4911EXP              0x00044911l
#define TYP_M2I4912EXP              0x00044912l
#define TYP_M2I4931EXP              0x00044931l
#define TYP_M2I4932EXP              0x00044932l
#define TYP_M2I4960EXP              0x00044960l
#define TYP_M2I4961EXP              0x00044961l
#define TYP_M2I4963EXP              0x00044963l
#define TYP_M2I4964EXP              0x00044964l

// ----- M2p.59xx-x4 -----
#define TYP_M2P59XX_X4              0x00095900l      // generic
#define TYP_M2P5911_X4              0x00095911l
#define TYP_M2P5912_X4              0x00095912l
#define TYP_M2P5913_X4              0x00095913l
#define TYP_M2P5916_X4              0x00095916l
#define TYP_M2P5920_X4              0x00095920l
#define TYP_M2P5921_X4              0x00095921l
#define TYP_M2P5922_X4              0x00095922l
#define TYP_M2P5923_X4              0x00095923l
#define TYP_M2P5926_X4              0x00095926l
#define TYP_M2P5930_X4              0x00095930l
#define TYP_M2P5931_X4              0x00095931l
#define TYP_M2P5932_X4              0x00095932l
#define TYP_M2P5933_X4              0x00095933l
#define TYP_M2P5936_X4              0x00095936l
#define TYP_M2P5940_X4              0x00095940l
#define TYP_M2P5941_X4              0x00095941l
#define TYP_M2P5942_X4              0x00095942l
#define TYP_M2P5943_X4              0x00095943l
#define TYP_M2P5946_X4              0x00095946l
#define TYP_M2P5960_X4              0x00095960l
#define TYP_M2P5961_X4              0x00095961l
#define TYP_M2P5962_X4              0x00095962l
#define TYP_M2P5966_X4              0x00095966l
#define TYP_M2P5968_X4              0x00095968l


// ----- MI.60xx, MC.60xx, MX.60xx -----
#define TYP_MI6010                  0x00006010l
#define TYP_MI6011                  0x00006011l
#define TYP_MI6012                  0x00006012l
#define TYP_MI6021                  0x00006021l
#define TYP_MI6022                  0x00006022l
#define TYP_MI6030                  0x00006030l
#define TYP_MI6031                  0x00006031l
#define TYP_MI6033                  0x00006033l
#define TYP_MI6034                  0x00006034l

#define TYP_M2I6010                 0x00036010l
#define TYP_M2I6011                 0x00036011l
#define TYP_M2I6012                 0x00036012l
#define TYP_M2I6021                 0x00036021l
#define TYP_M2I6022                 0x00036022l
#define TYP_M2I6030                 0x00036030l
#define TYP_M2I6031                 0x00036031l
#define TYP_M2I6033                 0x00036033l
#define TYP_M2I6034                 0x00036034l

#define TYP_M2I6010EXP              0x00046010l
#define TYP_M2I6011EXP              0x00046011l
#define TYP_M2I6012EXP              0x00046012l
#define TYP_M2I6021EXP              0x00046021l
#define TYP_M2I6022EXP              0x00046022l
#define TYP_M2I6030EXP              0x00046030l
#define TYP_M2I6031EXP              0x00046031l
#define TYP_M2I6033EXP              0x00046033l
#define TYP_M2I6034EXP              0x00046034l

#define TYP_MC6010                  0x00016010l
#define TYP_MC6011                  0x00016011l
#define TYP_MC6012                  0x00016012l
#define TYP_MC6021                  0x00016021l
#define TYP_MC6022                  0x00016022l
#define TYP_MC6030                  0x00016030l
#define TYP_MC6031                  0x00016031l
#define TYP_MC6033                  0x00016033l
#define TYP_MC6034                  0x00016034l

#define TYP_MX6010                  0x00026010l
#define TYP_MX6011                  0x00026011l
#define TYP_MX6021                  0x00026021l
#define TYP_MX6030                  0x00026030l
#define TYP_MX6033                  0x00026033l



// ----- MI.61xx, MC.61xx, MX.61xx -----
#define TYP_MI6105                  0x00006105l
#define TYP_MI6110                  0x00006110l
#define TYP_MI6111                  0x00006111l

#define TYP_M2I6105                 0x00036105l
#define TYP_M2I6110                 0x00036110l
#define TYP_M2I6111                 0x00036111l

#define TYP_M2I6105EXP              0x00046105l
#define TYP_M2I6110EXP              0x00046110l
#define TYP_M2I6111EXP              0x00046111l

#define TYP_MC6110                  0x00016110l
#define TYP_MC6111                  0x00016111l

#define TYP_MX6110                  0x00026110l

// ----- M2p.65xx-x4 -----
#define TYP_M2P65XX_X4              0x00096500l      // generic
#define TYP_M2P6522_X4              0x00096522l      // 4 ch @   40 MS/s (1x4) (low voltage)
#define TYP_M2P6523_X4              0x00096523l      // 8 ch @   40 MS/s (low voltage)
#define TYP_M2P6530_X4              0x00096530l      // 1 ch @   40 MS/s
#define TYP_M2P6531_X4              0x00096531l      // 2 ch @   40 MS/s
#define TYP_M2P6532_X4              0x00096532l      // 4 ch @   40 MS/s (1x4)
#define TYP_M2P6536_X4              0x00096536l      // 4 ch @   40 MS/s (2x2)
#define TYP_M2P6533_X4              0x00096533l      // 8 ch @   40 MS/s
#define TYP_M2P6540_X4              0x00096540l      // 1 ch @   40 MS/s (high voltage)
#define TYP_M2P6541_X4              0x00096541l      // 2 ch @   40 MS/s (high voltage)
#define TYP_M2P6546_X4              0x00096546l      // 4 ch @   40 MS/s (2x2) (high voltage)
#define TYP_M2P6560_X4              0x00096560l      // 1 ch @  125 MS/s
#define TYP_M2P6561_X4              0x00096561l      // 2 ch @  125 MS/s
#define TYP_M2P6562_X4              0x00096562l      // 4 ch @  125 MS/s (1x4)
#define TYP_M2P6566_X4              0x00096566l      // 4 ch @  125 MS/s (2x2)
#define TYP_M2P6568_X4              0x00096568l      // 8 ch @  125/80 MS/s
#define TYP_M2P6570_X4              0x00096570l      // 1 ch @  125 MS/s (high voltage)
#define TYP_M2P6571_X4              0x00096571l      // 2 ch @  125 MS/s (high voltage)
#define TYP_M2P6576_X4              0x00096576l      // 4 ch @  125 MS/s (2x2) (high voltage)

// ----- M4i.66xx-x8 (16 bit D/A) -----
// ----- M4i.66xx-x8 (generic) -----
#define TYP_M4I66XX_X8              0x00076600

#define TYP_M4I6620_X8              0x00076620      // 1 ch @  625 MS/s
#define TYP_M4I6621_X8              0x00076621      // 2 ch @  625 MS/s
#define TYP_M4I6622_X8              0x00076622      // 4 ch @  625 MS/s
#define TYP_M4I6630_X8              0x00076630      // 1 ch @ 1250 MS/s
#define TYP_M4I6631_X8              0x00076631      // 2 ch @ 1250 MS/s

// ----- M4x.66xx-x8 (16 bit D/A) -----
// ----- M4x.66xx-x8 (generic) -----
#define TYP_M4X66XX_X4              0x00086600

#define TYP_M4X6620_X4              0x00086620      // 1 ch @  625 MS/s
#define TYP_M4X6621_X4              0x00086621      // 2 ch @  625 MS/s
#define TYP_M4X6622_X4              0x00086622      // 4 ch @  625 MS/s
#define TYP_M4X6630_X4              0x00086630      // 1 ch @ 1250 MS/s
#define TYP_M4X6631_X4              0x00086631      // 2 ch @ 1250 MS/s

// ----- MI.70xx, MC.70xx, MX.70xx -----
#define TYP_MI7005                  0x00007005l
#define TYP_MI7010                  0x00007010l
#define TYP_MI7011                  0x00007011l
#define TYP_MI7020                  0x00007020l
#define TYP_MI7021                  0x00007021l

#define TYP_M2I7005                 0x00037005l
#define TYP_M2I7010                 0x00037010l
#define TYP_M2I7011                 0x00037011l
#define TYP_M2I7020                 0x00037020l
#define TYP_M2I7021                 0x00037021l

#define TYP_M2I7005EXP              0x00047005l
#define TYP_M2I7010EXP              0x00047010l
#define TYP_M2I7011EXP              0x00047011l
#define TYP_M2I7020EXP              0x00047020l
#define TYP_M2I7021EXP              0x00047021l

#define TYP_MC7005                  0x00017005l
#define TYP_MC7010                  0x00017010l
#define TYP_MC7011                  0x00017011l
#define TYP_MC7020                  0x00017020l
#define TYP_MC7021                  0x00017021l

#define TYP_MX7005                  0x00027005l
#define TYP_MX7010                  0x00027010l
#define TYP_MX7011                  0x00027011l



// ----- MI.72xx, MC.72xx, MX.72xx -----
#define TYP_MI7210                  0x00007210l
#define TYP_MI7211                  0x00007211l
#define TYP_MI7220                  0x00007220l
#define TYP_MI7221                  0x00007221l

#define TYP_M2I7210                 0x00037210l
#define TYP_M2I7211                 0x00037211l
#define TYP_M2I7220                 0x00037220l
#define TYP_M2I7221                 0x00037221l

#define TYP_M2I7210EXP              0x00047210l
#define TYP_M2I7211EXP              0x00047211l
#define TYP_M2I7220EXP              0x00047220l
#define TYP_M2I7221EXP              0x00047221l

#define TYP_MC7210                  0x00017210l
#define TYP_MC7211                  0x00017211l
#define TYP_MC7220                  0x00017220l
#define TYP_MC7221                  0x00017221l

#define TYP_MX7210                  0x00027210l
#define TYP_MX7220                  0x00027220l

// ----- M2p.75xx-x4 -----
#define TYP_M2P75XX_X4              0x00097500l      // generic
#define TYP_M2P7515_X4              0x00097515l

// ----- M4i.77xx-x8  -----
#define TYP_M4I77XX_X8              0x00077700 // generic
#define TYP_M4I7710_X8              0x00077710 // single-ended
#define TYP_M4I7720_X8              0x00077720 // single-ended
#define TYP_M4I7730_X8              0x00077730 // single-ended
#define TYP_M4I7725_X8              0x00077725 // differential
#define TYP_M4I7735_X8              0x00077735 // differential

// ----- M4x.77xx-x8  -----
#define TYP_M4X77XX_X4              0x00087700 // generic
#define TYP_M4X7710_X4              0x00087710 // single-ended
#define TYP_M4X7720_X4              0x00087720 // single-ended
#define TYP_M4X7730_X4              0x00087730 // single-ended
#define TYP_M4X7725_X4              0x00087725 // differential
#define TYP_M4X7735_X4              0x00087735 // differential

// ----- MX.90xx -----
#define TYP_MX9010                  0x00029010l



// ***********************************************************************
// software registers
// ***********************************************************************


// ***** PCI Features Bits (MI/MC/MX and prior cards) *********
#define PCIBIT_MULTI                0x00000001
#define PCIBIT_DIGITAL              0x00000002
#define PCIBIT_CH0DIGI              0x00000004
#define PCIBIT_EXTSAM               0x00000008
#define PCIBIT_3CHANNEL             0x00000010
#define PCIBIT_GATE                 0x00000020
#define PCIBIT_SLAVE                0x00000040
#define PCIBIT_MASTER               0x00000080
#define PCIBIT_DOUBLEMEM            0x00000100
#define PCIBIT_SYNC                 0x00000200
#define PCIBIT_TIMESTAMP            0x00000400
#define PCIBIT_STARHUB              0x00000800
#define PCIBIT_CA                   0x00001000
#define PCIBIT_XIO                  0x00002000
#define PCIBIT_AMPLIFIER            0x00004000
#define PCIBIT_DIFFMODE             0x00008000

#define PCIBIT_ELISA                0x10000000


// ***** PCI features starting with M2i card series *****
#define SPCM_FEAT_MULTI             0x00000001      // multiple recording
#define SPCM_FEAT_GATE              0x00000002      // gated sampling
#define SPCM_FEAT_DIGITAL           0x00000004      // additional synchronous digital inputs or outputs
#define SPCM_FEAT_TIMESTAMP         0x00000008      // timestamp
#define SPCM_FEAT_STARHUB5          0x00000020      // starhub for  5 cards installed (M2i + M2i-Exp)
#define SPCM_FEAT_STARHUB4          0x00000020      // starhub for  4 cards installed (M3i + M3i-Exp)
#define SPCM_FEAT_STARHUB6_EXTM     0x00000020      // starhub for  6 cards installed as card extension or piggy back (M2p)
#define SPCM_FEAT_STARHUB8_EXTM     0x00000020      // starhub for  8 cards installed as card extension or piggy back (M4i-Exp)
#define SPCM_FEAT_STARHUB16         0x00000040      // starhub for 16 cards installed (M2i, M2i-exp)
#define SPCM_FEAT_STARHUB16_EXTM    0x00000040      // starhub for 16 cards installed as card extension or piggy back (M2p)
#define SPCM_FEAT_STARHUB8          0x00000040      // starhub for  8 cards installed (M3i + M3i-Exp)
#define SPCM_FEAT_STARHUBXX_MASK    0x00000060      // mask to detect one of the above installed starhub
#define SPCM_FEAT_ABA               0x00000080      // ABA mode installed
#define SPCM_FEAT_BASEXIO           0x00000100      // extra I/O on base card installed
#define SPCM_FEAT_AMPLIFIER_10V     0x00000200      // external amplifier for 60/61
#define SPCM_FEAT_STARHUBSYSMASTER  0x00000400      // system starhub master installed
#define SPCM_FEAT_DIFFMODE          0x00000800      // Differential mode installed
#define SPCM_FEAT_SEQUENCE          0x00001000      // Sequence programming mode for generator cards
#define SPCM_FEAT_AMPMODULE_10V     0x00002000      // amplifier module for 60/61
#define SPCM_FEAT_STARHUBSYSSLAVE   0x00004000      // system starhub slave installed
#define SPCM_FEAT_NETBOX            0x00008000      // card is part of netbox
#define SPCM_FEAT_REMOTESERVER      0x00010000      // remote server can be used with this card
#define SPCM_FEAT_SCAPP             0x00020000      // SCAPP option (CUDA RDMA)
#define SPCM_FEAT_DIG16_SMB         0x00040000      // M2p: 16 additional digital inputs or outputs (via SMB connectors) 
#define SPCM_FEAT_DIG8_SMA          0x00040000      // M4i:  8 additional digital inputs or 6 additional outputs (via SMA connectors) 
#define SPCM_FEAT_DIG16_FX2         0x00080000      // M2p: 16 additional digital inputs or outputs (via FX2 connector)
#define SPCM_FEAT_DIGITALBWFILTER   0x00100000      // Digital BW filter is available
#define SPCM_FEAT_CUSTOMMOD_MASK    0xF0000000      // mask for custom modification code, meaning of code depends on type and customer


// ***** Extended Features starting with M4i *****
#define SPCM_FEAT_EXTFW_SEGSTAT     0x00000001        // segment (Multiple Recording, ABA) statistics like average, min/max
#define SPCM_FEAT_EXTFW_SEGAVERAGE  0x00000002        // average of multiple segments (Multiple Recording, ABA) 
#define SPCM_FEAT_EXTFW_BOXCAR      0x00000004      // boxcar averaging (high-res mode)


// ***** Error Request *************
#define ERRORTEXTLEN                200
#define SPC_LASTERRORTEXT           999996l
#define SPC_LASTERRORVALUE          999997l
#define SPC_LASTERRORREG            999998l
#define SPC_LASTERRORCODE           999999l     // Reading this reset the internal error-memory.

// ***** constants to use with the various _ACDC registers *****
#define COUPLING_DC 0
#define COUPLING_AC 1


// ***** Register and Command Structure
#define SPC_COMMAND                 0l
#define     SPC_RESET               0l
#define     SPC_SOFTRESET           1l
#define     SPC_WRITESETUP          2l
#define     SPC_START               10l
#define     SPC_STARTANDWAIT        11l
#define     SPC_FIFOSTART           12l
#define     SPC_FIFOWAIT            13l
#define     SPC_FIFOSTARTNOWAIT     14l
#define     SPC_FORCETRIGGER        16l
#define     SPC_STOP                20l
#define     SPC_FLUSHFIFOBUFFER     21l
#define     SPC_POWERDOWN           30l
#define     SPC_SYNCMASTER          100l
#define     SPC_SYNCTRIGGERMASTER   101l
#define     SPC_SYNCMASTERFIFO      102l
#define     SPC_SYNCSLAVE           110l
#define     SPC_SYNCTRIGGERSLAVE    111l
#define     SPC_SYNCSLAVEFIFO       112l
#define     SPC_NOSYNC              120l
#define     SPC_SYNCSTART           130l
#define     SPC_SYNCCALCMASTER      140l
#define     SPC_SYNCCALCMASTERFIFO  141l
#define     SPC_PXIDIVIDERRESET     150l
#define     SPC_RELAISON            200l
#define     SPC_RELAISOFF           210l
#define     SPC_ADJUSTSTART         300l
#define     SPC_FIFO_BUFREADY0      400l
#define     SPC_FIFO_BUFREADY1      401l
#define     SPC_FIFO_BUFREADY2      402l
#define     SPC_FIFO_BUFREADY3      403l
#define     SPC_FIFO_BUFREADY4      404l
#define     SPC_FIFO_BUFREADY5      405l
#define     SPC_FIFO_BUFREADY6      406l
#define     SPC_FIFO_BUFREADY7      407l
#define     SPC_FIFO_BUFREADY8      408l
#define     SPC_FIFO_BUFREADY9      409l
#define     SPC_FIFO_BUFREADY10     410l
#define     SPC_FIFO_BUFREADY11     411l
#define     SPC_FIFO_BUFREADY12     412l
#define     SPC_FIFO_BUFREADY13     413l
#define     SPC_FIFO_BUFREADY14     414l
#define     SPC_FIFO_BUFREADY15     415l
#define     SPC_FIFO_AUTOBUFSTART   500l
#define     SPC_FIFO_AUTOBUFEND     510l

#define SPC_STATUS                  10l
#define     SPC_RUN                 0l
#define     SPC_TRIGGER             10l
#define     SPC_READY               20l



// commands for M2 cards
#define SPC_M2CMD                   100l                // write a command
#define     M2CMD_CARD_RESET            0x00000001l     // hardware reset   
#define     M2CMD_CARD_WRITESETUP       0x00000002l     // write setup only
#define     M2CMD_CARD_START            0x00000004l     // start of card (including writesetup)
#define     M2CMD_CARD_ENABLETRIGGER    0x00000008l     // enable trigger engine
#define     M2CMD_CARD_FORCETRIGGER     0x00000010l     // force trigger
#define     M2CMD_CARD_DISABLETRIGGER   0x00000020l     // disable trigger engine again (multi or gate)
#define     M2CMD_CARD_STOP             0x00000040l     // stop run
#define     M2CMD_CARD_FLUSHFIFO        0x00000080l     // flush fifos to memory
#define     M2CMD_CARD_INVALIDATEDATA   0x00000100l     // current data in memory is invalidated, next data transfer start will wait until new data is available
#define     M2CMD_CARD_INTERNALRESET    0x00000200l     // INTERNAL reset command

#define     M2CMD_ALL_STOP              0x00440060l     // stops card and all running transfers

#define     M2CMD_CARD_WAITPREFULL      0x00001000l     // wait until pretrigger is full
#define     M2CMD_CARD_WAITTRIGGER      0x00002000l     // wait for trigger recognition
#define     M2CMD_CARD_WAITREADY        0x00004000l     // wait for card ready

#define     M2CMD_DATA_STARTDMA         0x00010000l     // start of DMA transfer for data
#define     M2CMD_DATA_WAITDMA          0x00020000l     // wait for end of data transfer / next block ready
#define     M2CMD_DATA_STOPDMA          0x00040000l     // abort the data transfer
#define     M2CMD_DATA_POLL             0x00080000l     // transfer data using single access and polling

#define     M2CMD_EXTRA_STARTDMA        0x00100000l     // start of DMA transfer for extra (ABA + timestamp) data
#define     M2CMD_EXTRA_WAITDMA         0x00200000l     // wait for end of extra (ABA + timestamp) data transfer / next block ready
#define     M2CMD_EXTRA_STOPDMA         0x00400000l     // abort the extra (ABA + timestamp) data transfer
#define     M2CMD_EXTRA_POLL            0x00800000l     // transfer data using single access and polling

#define     M2CMD_DATA_SGFLUSH          0x01000000l     // flush incomplete pages from sg list


// status for M2 cards (bitmask)
#define SPC_M2STATUS                110l                // read the current status
#define     M2STAT_NONE                 0x00000000l     // status empty
#define     M2STAT_CARD_PRETRIGGER      0x00000001l     // pretrigger area is full
#define     M2STAT_CARD_TRIGGER         0x00000002l     // trigger recognized
#define     M2STAT_CARD_READY           0x00000004l     // card is ready, run finished
#define     M2STAT_CARD_SEGMENT_PRETRG  0x00000008l     // since M4i: at muliple-recording: pretrigger area of a segment is full

#define     M2STAT_DATA_BLOCKREADY      0x00000100l     // next data block is available
#define     M2STAT_DATA_END             0x00000200l     // data transfer has ended
#define     M2STAT_DATA_OVERRUN         0x00000400l     // FIFO overrun (record) or underrun (replay)
#define     M2STAT_DATA_ERROR           0x00000800l     // internal error

#define     M2STAT_EXTRA_BLOCKREADY     0x00001000l     // next extra data (ABA and timestamp) block is available
#define     M2STAT_EXTRA_END            0x00002000l     // extra data (ABA and timestamp) transfer has ended
#define     M2STAT_EXTRA_OVERRUN        0x00004000l     // FIFO overrun
#define     M2STAT_EXTRA_ERROR          0x00008000l     // internal error

#define     M2STAT_TSCNT_OVERRUN        0x00010000l     // timestamp counter overrun

#define     M2STAT_INTERNALMASK         0xff000000l     // mask for internal status signals
#define     M2STAT_INTERNAL_SYSLOCK     0x02000000l



// buffer control registers for samples data
#define SPC_DATA_AVAIL_USER_LEN     200l                // number of bytes available for user (valid data if READ, free buffer if WRITE)
#define SPC_DATA_AVAIL_USER_POS     201l                // the current byte position where the available user data starts
#define SPC_DATA_AVAIL_CARD_LEN     202l                // number of bytes available for card (free buffer if READ, filled data if WRITE)
#define SPC_DATA_OUTBUFSIZE         209l                // output buffer size in bytes

// buffer control registers for extra data (ABA slow data, timestamps)
#define SPC_ABA_AVAIL_USER_LEN      210l                // number of bytes available for user (valid data if READ, free buffer if WRITE)
#define SPC_ABA_AVAIL_USER_POS      211l                // the current byte position where the available user data starts
#define SPC_ABA_AVAIL_CARD_LEN      212l                // number of bytes available for card (free buffer if READ, filled data if WRITE)

#define SPC_TS_AVAIL_USER_LEN       220l                // number of bytes available for user (valid data if READ, free buffer if WRITE)
#define SPC_TS_AVAIL_USER_POS       221l                // the current byte position where the available user data starts
#define SPC_TS_AVAIL_CARD_LEN       222l                // number of bytes available for card (free buffer if READ, filled data if WRITE)



// Installation
#define SPC_VERSION                 1000l
#define SPC_ISAADR                  1010l
#define SPC_INSTMEM                 1020l
#define SPC_INSTSAMPLERATE          1030l
#define SPC_BRDTYP                  1040l

// MI/MC/MX type information (internal use)
#define SPC_MIINST_MODULES          1100l
#define SPC_MIINST_CHPERMODULE      1110l
#define SPC_MIINST_BYTESPERSAMPLE   1120l
#define SPC_MIINST_BITSPERSAMPLE    1125l
#define SPC_MIINST_MAXADCVALUE      1126l
#define SPC_MIINST_MINADCLOCK       1130l
#define SPC_MIINST_MAXADCLOCK       1140l
#define SPC_MIINST_MINEXTCLOCK      1145l
#define SPC_MIINST_MAXEXTCLOCK      1146l
#define SPC_MIINST_MINSYNCCLOCK     1147l
#define SPC_MIINST_MINEXTREFCLOCK   1148l
#define SPC_MIINST_MAXEXTREFCLOCK   1149l
#define SPC_MIINST_QUARZ            1150l
#define SPC_MIINST_QUARZ2           1151l
#define SPC_MIINST_MINEXTCLOCK1     1152l
#define SPC_MIINST_FLAGS            1160l
#define SPC_MIINST_FIFOSUPPORT      1170l
#define SPC_MIINST_ISDEMOCARD       1175l

// Driver information
#define SPC_GETDRVVERSION           1200l
#define SPC_GETKERNELVERSION        1210l
#define SPC_GETDRVTYPE              1220l
#define     DRVTYP_DOS              0l
#define     DRVTYP_LINUX32          1l
#define     DRVTYP_VXD              2l
#define     DRVTYP_NTLEGACY         3l
#define     DRVTYP_WDM32            4l
#define     DRVTYP_WDM64            5l
#define     DRVTYP_WOW64            6l
#define     DRVTYP_LINUX64          7l
#define     DRVTYP_QNX32            8l
#define     DRVTYP_QNX64            9l
#define SPC_GETCOMPATIBILITYVERSION 1230l
#define SPC_GETMINDRVVERSION        1240l

// PCI, CompactPCI and PXI Installation Information
#define SPC_PCITYP                  2000l

// ***** available card function types *****
#define SPC_FNCTYPE                 2001l
#define     SPCM_TYPE_AI            0x01
#define     SPCM_TYPE_AO            0x02
#define     SPCM_TYPE_DI            0x04
#define     SPCM_TYPE_DO            0x08
#define     SPCM_TYPE_DIO           0x10

#define SPC_PCIVERSION              2010l
#define SPC_PCIEXTVERSION           2011l
#define SPC_PCIMODULEVERSION        2012l
#define SPC_PCIMODULEBVERSION       2013l
#define SPC_BASEPCBVERSION          2014l
#define SPC_MODULEPCBVERSION        2015l
#define SPC_MODULEAPCBVERSION       2015l
#define SPC_MODULEBPCBVERSION       2016l
#define SPC_EXTPCBVERSION           2017l
#define SPC_PCIDIGVERSION           2018l
#define SPC_DIGPCBVERSION           2019l
#define SPC_PCIDATE                 2020l
#define SPC_CALIBDATE               2025l
#define SPC_CALIBDATEONBOARD        2026l
#define SPC_PCISERIALNR             2030l
#define SPC_PCISERIALNO             2030l
#define SPC_PCIHWBUSNO              2040l
#define SPC_PCIHWDEVNO              2041l
#define SPC_PCIHWFNCNO              2042l
#define SPC_PCIHWSLOTNO             2043l
#define SPC_PCIEXPGENERATION        2050l
#define SPC_PCIEXPLANES             2051l
#define SPC_PCIEXPPAYLOAD           2052l
#define SPC_PCIEXPREADREQUESTSIZE   2053l
#define SPC_PCIEXPREADCOMPLBOUNDARY 2054l
#define SPC_PXIHWSLOTNO             2055l
#define SPC_PCISAMPLERATE           2100l
#define SPC_PCIMEMSIZE              2110l
#define SPC_PCIFEATURES             2120l
#define SPC_PCIEXTFEATURES          2121l
#define SPC_PCIINFOADR              2200l
#define SPC_PCIINTERRUPT            2300l
#define SPC_PCIBASEADR0             2400l
#define SPC_PCIBASEADR1             2401l
#define SPC_PCIREGION0              2410l
#define SPC_PCIREGION1              2411l
#define SPC_READTRGLVLCOUNT         2500l
#define SPC_READIRCOUNT             3000l
#define SPC_READUNIPOLAR0           3010l
#define SPC_READUNIPOLAR1           3020l
#define SPC_READUNIPOLAR2           3030l
#define SPC_READUNIPOLAR3           3040l
#define SPC_READMAXOFFSET           3100l

#define SPC_READAIFEATURES          3101l
#define     SPCM_AI_TERM            0x00000001  // input termination available
#define     SPCM_AI_SE              0x00000002  // single-ended mode available
#define     SPCM_AI_DIFF            0x00000004  // differential mode available
#define     SPCM_AI_OFFSPERCENT     0x00000008  // offset programming is done in percent of input range
#define     SPCM_AI_OFFSMV          0x00000010  // offset programming is done in mV absolut
#define     SPCM_AI_OVERRANGEDETECT 0x00000020  // overrange detection is programmable
#define     SPCM_AI_DCCOUPLING      0x00000040  // DC coupling available
#define     SPCM_AI_ACCOUPLING      0x00000080  // AC coupling available
#define     SPCM_AI_LOWPASS         0x00000100  // selectable low pass
#define     SPCM_AI_ACDC_OFFS_COMP  0x00000200  // AC/DC offset compensation
#define     SPCM_AI_DIFFMUX         0x00000400  // differential mode (two channels combined to one) available
#define     SPCM_AI_GLOBALLOWPASS   0x00000800  // globally selectable low pass (all channels same setting)
#define     SPCM_AI_AUTOCALOFFS     0x00001000  // automatic offset calibration in hardware
#define     SPCM_AI_AUTOCALGAIN     0x00002000  // automatic gain calibration in hardware
#define     SPCM_AI_AUTOCALOFFSNOIN 0x00004000  // automatic offset calibration with open inputs
#define     SPCM_AI_HIGHIMP         0x00008000  // high impedance available
#define     SPCM_AI_LOWIMP          0x00010000  // low impedance available (50 ohm)
#define     SPCM_AI_DIGITALLOWPASS  0x00020000  // selectable digital low pass filter
#define     SPCM_AI_INDIVPULSEWIDTH 0x00100000  // individual pulsewidth per channel available

#define SPC_READAOFEATURES          3102l
#define     SPCM_AO_SE              0x00000002  // single-ended mode available
#define     SPCM_AO_DIFF            0x00000004  // differential mode available
#define     SPCM_AO_PROGFILTER      0x00000008  // programmable filters available
#define     SPCM_AO_PROGOFFSET      0x00000010  // programmable offset available
#define     SPCM_AO_PROGGAIN        0x00000020  // programmable gain available
#define     SPCM_AO_PROGSTOPLEVEL   0x00000040  // programmable stop level available
#define     SPCM_AO_DOUBLEOUT       0x00000080  // double out mode available
#define     SPCM_AO_ENABLEOUT       0x00000100  // outputs can be disabled/enabled

#define SPC_READDIFEATURES          3103l
#define     SPCM_DI_TERM            0x00000001  // input termination available
#define     SPCM_DI_SE              0x00000002  // single-ended mode available
#define     SPCM_DI_DIFF            0x00000004  // differential mode available
#define     SPCM_DI_PROGTHRESHOLD   0x00000008  // programmable threshold available
#define     SPCM_DI_HIGHIMP         0x00000010  // high impedance available
#define     SPCM_DI_LOWIMP          0x00000020  // low impedance available
#define     SPCM_DI_INDIVPULSEWIDTH 0x00100000  // individual pulsewidth per channel available
#define     SPCM_DI_IOCHANNEL       0x00200000  // connected with DO channel

#define SPC_READDOFEATURES          3104l
#define     SPCM_DO_SE              0x00000002  // single-ended mode available
#define     SPCM_DO_DIFF            0x00000004  // differential mode available
#define     SPCM_DO_PROGSTOPLEVEL   0x00000008  // programmable stop level available
#define     SPCM_DO_PROGOUTLEVELS   0x00000010  // programmable output levels (low + high) available
#define     SPCM_DO_ENABLEMASK      0x00000020  // individual enable mask for each output channel
#define     SPCM_DO_IOCHANNEL       0x00200000  // connected with DI channel

#define SPC_READCHGROUPING          3110l
#define SPC_READAIPATHCOUNT         3120l       // number of available analog input paths
#define SPC_READAIPATH              3121l       // the current path for which all the settings are read

#define SPCM_CUSTOMMOD              3130l
#define     SPCM_CUSTOMMOD_BASE_MASK    0x000000FF
#define     SPCM_CUSTOMMOD_MODULE_MASK  0x0000FF00
#define     SPCM_CUSTOMMOD_STARHUB_MASK 0x00FF0000

#define SPC_READRANGECH0_0          3200l
#define SPC_READRANGECH0_1          3201l
#define SPC_READRANGECH0_2          3202l
#define SPC_READRANGECH0_3          3203l
#define SPC_READRANGECH0_4          3204l
#define SPC_READRANGECH0_5          3205l
#define SPC_READRANGECH0_6          3206l
#define SPC_READRANGECH0_7          3207l
#define SPC_READRANGECH0_8          3208l
#define SPC_READRANGECH0_9          3209l
#define SPC_READRANGECH1_0          3300l
#define SPC_READRANGECH1_1          3301l
#define SPC_READRANGECH1_2          3302l
#define SPC_READRANGECH1_3          3303l
#define SPC_READRANGECH1_4          3304l
#define SPC_READRANGECH1_5          3305l
#define SPC_READRANGECH1_6          3306l
#define SPC_READRANGECH1_7          3307l
#define SPC_READRANGECH1_8          3308l
#define SPC_READRANGECH1_9          3309l
#define SPC_READRANGECH2_0          3400l
#define SPC_READRANGECH2_1          3401l
#define SPC_READRANGECH2_2          3402l
#define SPC_READRANGECH2_3          3403l
#define SPC_READRANGECH3_0          3500l
#define SPC_READRANGECH3_1          3501l
#define SPC_READRANGECH3_2          3502l
#define SPC_READRANGECH3_3          3503l

#define SPC_READRANGEMIN0           4000l
#define SPC_READRANGEMIN99          4099l
#define SPC_READRANGEMAX0           4100l
#define SPC_READRANGEMAX99          4199l
#define SPC_READOFFSMIN0            4200l
#define SPC_READOFFSMIN99           4299l
#define SPC_READOFFSMAX0            4300l
#define SPC_READOFFSMAX99           4399l
#define SPC_PCICOUNTER              9000l
#define SPC_BUFFERPOS               9010l

#define SPC_READAOGAINMIN           9100l
#define SPC_READAOGAINMAX           9110l
#define SPC_READAOOFFSETMIN         9120l
#define SPC_READAOOFFSETMAX         9130l

#define SPC_CARDMODE                9500l       // card modes as listed below
#define SPC_AVAILCARDMODES          9501l       // list with available card modes

// card modes
#define     SPC_REC_STD_SINGLE          0x00000001  // singleshot recording to memory
#define     SPC_REC_STD_MULTI           0x00000002  // multiple records to memory on each trigger event
#define     SPC_REC_STD_GATE            0x00000004  // gated recording to memory on gate signal
#define     SPC_REC_STD_ABA             0x00000008  // ABA: A slowly to extra FIFO, B to memory on each trigger event 
#define     SPC_REC_STD_SEGSTATS        0x00010000  // segment information stored on each trigger segment -> stored in on-board memory
#define     SPC_REC_STD_AVERAGE         0x00020000  // multiple records summed to average memory on each trigger event -> stored in on-board memory
#define     SPC_REC_STD_AVERAGE_16BIT   0x00080000  // multiple records summed to average memory on each trigger event -> stored in on-board memory
#define     SPC_REC_STD_BOXCAR          0x00800000  // boxcar averaging

#define     SPC_REC_FIFO_SINGLE         0x00000010  // singleshot to FIFO on trigger event
#define     SPC_REC_FIFO_MULTI          0x00000020  // multiple records to FIFO on each trigger event
#define     SPC_REC_FIFO_GATE           0x00000040  // gated sampling to FIFO on gate signal
#define     SPC_REC_FIFO_ABA            0x00000080  // ABA: A slowly to extra FIFO, B to FIFO on each trigger event
#define     SPC_REC_FIFO_SEGSTATS       0x00100000  // segment information stored on each trigger segment -> streamed to host
#define     SPC_REC_FIFO_AVERAGE        0x00200000  // multiple records summed to average memory on each trigger event -> streamed to host
#define     SPC_REC_FIFO_AVERAGE_16BIT  0x00400000  // multiple records summed to average memory on each trigger event -> streamed to host
#define     SPC_REC_FIFO_BOXCAR         0x01000000  // boxcar averaging FIFO mode
#define     SPC_REC_FIFO_SINGLE_MONITOR 0x02000000  // like SPC_REC_FIFO_SINGLE but with additional slow A data stream for monitoring

#define     SPC_REP_STD_SINGLE          0x00000100  // single replay from memory on trigger event 
#define     SPC_REP_STD_MULTI           0x00000200  // multiple replay from memory on each trigger event
#define     SPC_REP_STD_GATE            0x00000400  // gated replay from memory on gate signal

#define     SPC_REP_FIFO_SINGLE         0x00000800  // single replay from FIFO on trigger event
#define     SPC_REP_FIFO_MULTI          0x00001000  // multiple replay from FIFO on each trigger event
#define     SPC_REP_FIFO_GATE           0x00002000  // gated replay from FIFO on gate signal

#define     SPC_REP_STD_CONTINUOUS      0x00004000  // continuous replay started by one trigger event
#define     SPC_REP_STD_SINGLERESTART   0x00008000  // single replays on every detected trigger event
#define     SPC_REP_STD_SEQUENCE        0x00040000  // sequence mode replay

// Waveforms for demo cards
#define SPC_DEMOWAVEFORM            9600l
#define SPC_AVAILDEMOWAVEFORMS      9601l
#define     SPCM_DEMOWAVEFORM_SINE      0x00000001
#define     SPCM_DEMOWAVEFORM_RECT      0x00000002
#define     SPCM_DEMOWAVEFORM_TRIANGLE  0x00000004


// Memory
#define SPC_MEMSIZE                 10000l
#define SPC_SEGMENTSIZE             10010l
#define SPC_LOOPS                   10020l
#define SPC_PRETRIGGER              10030l
#define SPC_ABADIVIDER              10040l
#define SPC_AVERAGES                10050l
#define SPC_BOX_AVERAGES            10060l
#define SPC_SEGSPLIT_START          10070l
#define SPC_SEGSPLIT_PAUSE          10071l
#define SPC_POSTTRIGGER             10100l
#define SPC_STARTOFFSET             10200l

// Memory info (depends on mode and channelenable)
#define SPC_AVAILMEMSIZE_MIN        10201l
#define SPC_AVAILMEMSIZE_MAX        10202l
#define SPC_AVAILMEMSIZE_STEP       10203l
#define SPC_AVAILPOSTTRIGGER_MIN    10204l
#define SPC_AVAILPOSTTRIGGER_MAX    10205l
#define SPC_AVAILPOSTTRIGGER_STEP   10206l

#define SPC_AVAILABADIVIDER_MIN     10207l
#define SPC_AVAILABADIVIDER_MAX     10208l
#define SPC_AVAILABADIVIDER_STEP    10209l

#define SPC_AVAILLOOPS_MIN          10210l
#define SPC_AVAILLOOPS_MAX          10211l
#define SPC_AVAILLOOPS_STEP         10212l

#define SPC_AVAILAVERAGES_MIN       10220l
#define SPC_AVAILAVERAGES_MAX       10221l
#define SPC_AVAILAVERAGES_STEP      10222l

#define SPC_AVAILAVRGSEGSIZE_MIN    10223l
#define SPC_AVAILAVRGSEGSIZE_MAX    10224l
#define SPC_AVAILAVRGSEGSIZE_STEP   10225l

#define SPC_AVAILAVERAGES16BIT_MIN     10226l
#define SPC_AVAILAVERAGES16BIT_MAX     10227l
#define SPC_AVAILAVERAGES16BIT_STEP    10228l

#define SPC_AVAILAVRG16BITSEGSIZE_MIN  10229l
#define SPC_AVAILAVRG16BITSEGSIZE_MAX  10230l
#define SPC_AVAILAVRG16BITSEGSIZE_STEP 10231l

#define SPC_AVAILBOXCARAVERAGES_MIN         10232l
#define SPC_AVAILBOXCARAVERAGES_MAX         10233l
#define SPC_AVAILBOXCARAVERAGES_STEPFACTOR  10234l


// Channels
#define SPC_CHENABLE                11000l
#define SPC_CHCOUNT                 11001l
#define SPC_CHMODACOUNT             11100l
#define SPC_CHMODBCOUNT             11101l


// ----- channel enable flags for A/D and D/A boards (MI/MC/MX series) -----
//       and all cards on M2i series
#define     CHANNEL0                0x00000001
#define     CHANNEL1                0x00000002
#define     CHANNEL2                0x00000004
#define     CHANNEL3                0x00000008
#define     CHANNEL4                0x00000010
#define     CHANNEL5                0x00000020
#define     CHANNEL6                0x00000040
#define     CHANNEL7                0x00000080
#define     CHANNEL8                0x00000100
#define     CHANNEL9                0x00000200
#define     CHANNEL10               0x00000400
#define     CHANNEL11               0x00000800
#define     CHANNEL12               0x00001000
#define     CHANNEL13               0x00002000
#define     CHANNEL14               0x00004000
#define     CHANNEL15               0x00008000
#define     CHANNEL16               0x00010000
#define     CHANNEL17               0x00020000
#define     CHANNEL18               0x00040000
#define     CHANNEL19               0x00080000
#define     CHANNEL20               0x00100000
#define     CHANNEL21               0x00200000
#define     CHANNEL22               0x00400000
#define     CHANNEL23               0x00800000
#define     CHANNEL24               0x01000000
#define     CHANNEL25               0x02000000
#define     CHANNEL26               0x04000000
#define     CHANNEL27               0x08000000
#define     CHANNEL28               0x10000000
#define     CHANNEL29               0x20000000
#define     CHANNEL30               0x40000000
#define     CHANNEL31               0x80000000
// CHANNEL32 up to CHANNEL63 are placed in the upper 32 bit of a 64 bit word (M2i only)


// ----- old digital i/o settings for 16 bit implementation (MI/MC/MX series)  -----
#define     CH0_8BITMODE            65536l  // for MI.70xx only
#define     CH0_16BIT               1l
#define     CH0_32BIT               3l
#define     CH1_16BIT               4l
#define     CH1_32BIT               12l

// ----- new digital i/o settings for 8 bit implementation (MI/MC/MX series) -----
#define     MOD0_8BIT               1l
#define     MOD0_16BIT              3l
#define     MOD0_32BIT              15l
#define     MOD1_8BIT               16l
#define     MOD1_16BIT              48l
#define     MOD1_32BIT              240l

#define SPC_CHROUTE0                11010l
#define SPC_CHROUTE1                11020l

#define SPC_BITENABLE               11030l



// ----- Clock Settings -----
#define SPC_SAMPLERATE              20000l
#define SPC_SYNCCLOCK               20005l
#define SPC_SAMPLERATE2             20010l
#define SPC_SR2                     20020l
#define SPC_PLL_ENABLE              20030l
#define SPC_PLL_ISLOCKED            20031l
#define SPC_CLOCKDIV                20040l
#define SPC_INTCLOCKDIV             20041l
#define SPC_PXICLOCKDIV             20042l
#define SPC_PLL_R                   20060l
#define SPC_PLL_F                   20061l
#define SPC_PLL_S                   20062l
#define SPC_PLL_DIV                 20063l
#define SPC_PXI_CLK_OUT             20090l
#define SPC_EXTERNALCLOCK           20100l
#define SPC_EXTERNOUT               20110l
#define SPC_CLOCKOUT                20110l
#define SPC_CLOCKOUTFREQUENCY       20111l
#define SPC_CLOCK50OHM              20120l
#define SPC_CLOCK110OHM             20120l
#define SPC_CLOCK75OHM              20120l
#define SPC_STROBE75OHM             20121l
#define SPC_EXTERNRANGE             20130l
#define SPC_EXTRANGESHDIRECT        20131l
#define     EXRANGE_NONE            0l
#define     EXRANGE_NOPLL           1l
#define     EXRANGE_SINGLE          2l
#define     EXRANGE_BURST_S         4l
#define     EXRANGE_BURST_M         8l
#define     EXRANGE_BURST_L         16l
#define     EXRANGE_BURST_XL        32l
#define     EXRANGE_LOW             64l
#define     EXRANGE_HIGH            128l
#define     EXRANGE_LOW_DPS         256l            // digital phase synchronization
#define SPC_REFERENCECLOCK          20140l
#define     REFCLOCK_PXI            -1l

// ----- new clock registers starting with M2i cards -----
#define SPC_CLOCKMODE               20200l      // clock mode as listed below
#define SPC_AVAILCLOCKMODES         20201l      // returns all available clock modes
#define     SPC_CM_INTPLL           0x00000001      // use internal PLL
#define     SPC_CM_QUARTZ1          0x00000002      // use plain quartz1 (with divider)
#define     SPC_CM_QUARTZ2          0x00000004      // use plain quartz2 (with divider)
#define     SPC_CM_EXTERNAL         0x00000008      // use external clock directly
#define     SPC_CM_EXTERNAL0        0x00000008      // use external clock0 directly (identical value to SPC_CM_EXTERNAL)
#define     SPC_CM_EXTDIVIDER       0x00000010      // use external clock with programmed divider
#define     SPC_CM_EXTREFCLOCK      0x00000020      // external reference clock fed in (defined with SPC_REFERENCECLOCK)
#define     SPC_CM_PXIREFCLOCK      0x00000040      // PXI reference clock
#define     SPC_CM_SHDIRECT         0x00000080      // Star-hub direct clock (not synchronised)
#define     SPC_CM_QUARTZ2_DIRSYNC  0x00000100      // use plain quartz2 (with divider) and put the Q2 clock on the star-hub module
#define     SPC_CM_QUARTZ1_DIRSYNC  0x00000200      // use plain quartz1 (with divider) and put the Q1 clock on the star-hub module
#define     SPC_CM_EXTERNAL1        0x00000400      // use external clock1 directly
// ----- internal use only! -----
#define     SPC_CM_SYNCINT          0x01000000
#define     SPC_CM_SYNCEXT          0x02000000

#define SPC_CLOCK_READFEATURES      20205l
#define SPC_CLOCK_READFEATURES0     20205l
#define SPC_CLOCK_READFEATURES1     20206l
#define     SPCM_CKFEAT_TERM            0x00000001
#define     SPCM_CKFEAT_HIGHIMP         0x00000002
#define     SPCM_CKFEAT_DCCOUPLING      0x00000004
#define     SPCM_CKFEAT_ACCOUPLING      0x00000008
#define     SPCM_CKFEAT_SE              0x00000010
#define     SPCM_CKFEAT_DIFF            0x00000020
#define     SPCM_CKFEAT_PROGEDGE        0x00000040
#define     SPCM_CKFEAT_LEVELPROG       0x00000100
#define     SPCM_CKFEAT_PROGTHRESHOLD   0x00000200
#define     SPCM_CKFEAT_PROGDELAY       0x00000400

#define SPC_BURSTSYSCLOCKMODE       20210l
#define SPC_SYNCMASTERSYSCLOCKMODE  20211l
#define SPC_CLOCK_SETUP_CHANGED     20212l

// clock delay if available
#define SPC_CLOCK_AVAILDELAY_MIN    20220l
#define SPC_CLOCK_AVAILDELAY_MAX    20221l
#define SPC_CLOCK_AVAILDELAY_STEP   20222l
#define SPC_CLOCK_DELAY             20223l

// clock edges
#define SPC_AVAILCLOCKEDGES         20224l
#define     SPCM_EDGE_FALLING       0x00000001 // Originally SPCM_RISING_EDGE  : name and value of constant intentionally changed with driver versions greater than V5.24. See hardware manual for details.
#define     SPCM_EDGE_RISING        0x00000002 // Originally SPCM_FALLING_EDGE : name and value of constant intentionally changed with driver versions greater than V5.24. See hardware manual for details.
#define     SPCM_BOTH_EDGES         0x00000004
#define     SPCM_EDGES_BOTH         0x00000004 //Just added for good measure to match naming scheme of above SPCM_EDGE_FALLING and SPCM_EDGE_RISING constants.
#define SPC_CLOCK_EDGE              20225l

// mux definitions for channel routing
#define SPC_CHANNELMUXINFO          20300l
#define     SPCM_MUX_NONE            0x00000000  // nothing is interlaced
#define     SPCM_MUX_MUXONMOD        0x00000001  // data on module is multiplexed, only one channel can have full speed
#define     SPCM_MUX_INVERTCLKONMOD  0x00000002  // two channels on one module run with inverted clock
#define     SPCM_MUX_DLY             0x00000003  // delay cable between modules, one channel can have full interlace speed
#define     SPCM_MUX_DLYANDMUXONMOD  0x00000004  // delay cable between modules and multplexing on module
#define     SPCM_MUX_MUXBETWEENMODS  0x00000005  // multiplexed between modules (fastest sampling rate only with one module)
#define     SPCM_MUX_MUXONMOD2CH     0x00000006  // data on module is multiplexed, only two channel can have full speed
#define     SPCM_MUX_MAX4CH          0x00000007  // only four channels can have full speed, independent of distribution on modules


// ----- In/Out Range -----
#define SPC_OFFS0                   30000l
#define SPC_AMP0                    30010l
#define SPC_ACDC0                   30020l
#define SPC_ACDC_OFFS_COMPENSATION0 30021l
#define SPC_50OHM0                  30030l
#define SPC_DIFF0                   30040l
#define SPC_DOUBLEOUT0              30041l
#define SPC_DIGITAL0                30050l
#define SPC_110OHM0                 30060l
#define SPC_110OHM0L                30060l
#define SPC_75OHM0                  30060l
#define SPC_INOUT0                  30070l
#define SPC_FILTER0                 30080l
#define SPC_BANKSWITCH0             30081l
#define SPC_PATH0                   30090l
#define SPC_ENABLEOUT0              30091l

#define SPC_OFFS1                   30100l
#define SPC_AMP1                    30110l
#define SPC_ACDC1                   30120l
#define SPC_ACDC_OFFS_COMPENSATION1 30121l
#define SPC_50OHM1                  30130l
#define SPC_DIFF1                   30140l
#define SPC_DOUBLEOUT1              30141l
#define SPC_DIGITAL1                30150l
#define SPC_110OHM1                 30160l
#define SPC_110OHM0H                30160l
#define SPC_75OHM1                  30160l
#define SPC_INOUT1                  30170l
#define SPC_FILTER1                 30180l
#define SPC_BANKSWITCH1             30181l
#define SPC_PATH1                   30190l
#define SPC_ENABLEOUT1              30191l

#define SPC_OFFS2                   30200l
#define SPC_AMP2                    30210l
#define SPC_ACDC2                   30220l
#define SPC_ACDC_OFFS_COMPENSATION2 30221l
#define SPC_50OHM2                  30230l
#define SPC_DIFF2                   30240l
#define SPC_DOUBLEOUT2              30241l
#define SPC_110OHM2                 30260l
#define SPC_110OHM1L                30260l
#define SPC_75OHM2                  30260l
#define SPC_INOUT2                  30270l
#define SPC_FILTER2                 30280l
#define SPC_BANKSWITCH2             30281l
#define SPC_PATH2                   30290l
#define SPC_ENABLEOUT2              30291l

#define SPC_OFFS3                   30300l
#define SPC_AMP3                    30310l
#define SPC_ACDC3                   30320l
#define SPC_ACDC_OFFS_COMPENSATION3 30321l
#define SPC_50OHM3                  30330l
#define SPC_DIFF3                   30340l
#define SPC_DOUBLEOUT3              30341l
#define SPC_110OHM3                 30360l
#define SPC_110OHM1H                30360l
#define SPC_75OHM3                  30360l
#define SPC_INOUT3                  30370l
#define SPC_FILTER3                 30380l
#define SPC_BANKSWITCH3             30381l
#define SPC_PATH3                   30390l
#define SPC_ENABLEOUT3              30391l

#define SPC_OFFS4                   30400l
#define SPC_AMP4                    30410l
#define SPC_ACDC4                   30420l
#define SPC_50OHM4                  30430l
#define SPC_DIFF4                   30440l
#define SPC_DOUBLEOUT4              30441l
#define SPC_FILTER4                 30480l
#define SPC_ENABLEOUT4              30491l
#define SPC_PATH4                   30490l

#define SPC_OFFS5                   30500l
#define SPC_AMP5                    30510l
#define SPC_ACDC5                   30520l
#define SPC_50OHM5                  30530l
#define SPC_DIFF5                   30540l
#define SPC_DOUBLEOUT5              30541l
#define SPC_FILTER5                 30580l
#define SPC_ENABLEOUT5              30591l
#define SPC_PATH5                   30590l

#define SPC_OFFS6                   30600l
#define SPC_AMP6                    30610l
#define SPC_ACDC6                   30620l
#define SPC_50OHM6                  30630l
#define SPC_DIFF6                   30640l
#define SPC_DOUBLEOUT6              30641l
#define SPC_FILTER6                 30680l
#define SPC_ENABLEOUT6              30691l
#define SPC_PATH6                   30690l

#define SPC_OFFS7                   30700l
#define SPC_AMP7                    30710l
#define SPC_ACDC7                   30720l
#define SPC_50OHM7                  30730l
#define SPC_DIFF7                   30740l
#define SPC_DOUBLEOUT7              30741l
#define SPC_FILTER7                 30780l
#define SPC_ENABLEOUT7              30791l
#define SPC_PATH7                   30790l

#define SPC_OFFS8                   30800l
#define SPC_AMP8                    30810l
#define SPC_ACDC8                   30820l
#define SPC_50OHM8                  30830l
#define SPC_DIFF8                   30840l
#define SPC_PATH8                   30890l

#define SPC_OFFS9                   30900l
#define SPC_AMP9                    30910l
#define SPC_ACDC9                   30920l
#define SPC_50OHM9                  30930l
#define SPC_DIFF9                   30940l
#define SPC_PATH9                   30990l

#define SPC_OFFS10                  31000l
#define SPC_AMP10                   31010l
#define SPC_ACDC10                  31020l
#define SPC_50OHM10                 31030l
#define SPC_DIFF10                  31040l
#define SPC_PATH10                  31090l

#define SPC_OFFS11                  31100l
#define SPC_AMP11                   31110l
#define SPC_ACDC11                  31120l
#define SPC_50OHM11                 31130l
#define SPC_DIFF11                  31140l
#define SPC_PATH11                  31190l

#define SPC_OFFS12                  31200l
#define SPC_AMP12                   31210l
#define SPC_ACDC12                  31220l
#define SPC_50OHM12                 31230l
#define SPC_DIFF12                  31240l
#define SPC_PATH12                  31290l

#define SPC_OFFS13                  31300l
#define SPC_AMP13                   31310l
#define SPC_ACDC13                  31320l
#define SPC_50OHM13                 31330l
#define SPC_DIFF13                  31340l
#define SPC_PATH13                  31390l

#define SPC_OFFS14                  31400l
#define SPC_AMP14                   31410l
#define SPC_ACDC14                  31420l
#define SPC_50OHM14                 31430l
#define SPC_DIFF14                  31440l
#define SPC_PATH14                  31490l

#define SPC_OFFS15                  31500l
#define SPC_AMP15                   31510l
#define SPC_ACDC15                  31520l
#define SPC_50OHM15                 31530l
#define SPC_DIFF15                  31540l
#define SPC_PATH15                  31590l

#define SPC_110OHMTRIGGER           30400l
#define SPC_110OHMCLOCK             30410l


#define   AMP_BI200                 200l
#define   AMP_BI500                 500l
#define   AMP_BI1000                1000l
#define   AMP_BI2000                2000l
#define   AMP_BI2500                2500l
#define   AMP_BI4000                4000l
#define   AMP_BI5000                5000l
#define   AMP_BI10000               10000l
#define   AMP_UNI400                100400l
#define   AMP_UNI1000               101000l
#define   AMP_UNI2000               102000l


// ----- Trigger Settings -----
#define SPC_TRIGGERMODE             40000l
#define SPC_TRIG_OUTPUT             40100l
#define SPC_TRIGGEROUT              40100l
#define SPC_TRIG_TERM               40110l
#define SPC_TRIG_TERM0              40110l
#define SPC_TRIGGER50OHM            40110l
#define SPC_TRIGGER110OHM0          40110l
#define SPC_TRIGGER75OHM0           40110l
#define SPC_TRIG_TERM1              40111l
#define SPC_TRIGGER110OHM1          40111l
#define SPC_TRIG_EXT0_ACDC          40120l
#define SPC_TRIG_EXT1_ACDC          40121l
#define SPC_TRIG_EXT2_ACDC          40122l

#define SPC_TRIGGERMODE0            40200l
#define SPC_TRIGGERMODE1            40201l
#define SPC_TRIGGERMODE2            40202l
#define SPC_TRIGGERMODE3            40203l
#define SPC_TRIGGERMODE4            40204l
#define SPC_TRIGGERMODE5            40205l
#define SPC_TRIGGERMODE6            40206l
#define SPC_TRIGGERMODE7            40207l
#define SPC_TRIGGERMODE8            40208l
#define SPC_TRIGGERMODE9            40209l
#define SPC_TRIGGERMODE10           40210l
#define SPC_TRIGGERMODE11           40211l
#define SPC_TRIGGERMODE12           40212l
#define SPC_TRIGGERMODE13           40213l
#define SPC_TRIGGERMODE14           40214l
#define SPC_TRIGGERMODE15           40215l

#define     TM_SOFTWARE             0l
#define     TM_NOTRIGGER            10l
#define     TM_CHXPOS               10000l
#define     TM_CHXPOS_LP            10001l
#define     TM_CHXPOS_SP            10002l
#define     TM_CHXPOS_GS            10003l
#define     TM_CHXPOS_SS            10004l
#define     TM_CHXNEG               10010l
#define     TM_CHXNEG_LP            10011l
#define     TM_CHXNEG_SP            10012l
#define     TM_CHXNEG_GS            10013l
#define     TM_CHXNEG_SS            10014l
#define     TM_CHXOFF               10020l
#define     TM_CHXBOTH              10030l
#define     TM_CHXWINENTER          10040l
#define     TM_CHXWINENTER_LP       10041l
#define     TM_CHXWINENTER_SP       10042l
#define     TM_CHXWINLEAVE          10050l
#define     TM_CHXWINLEAVE_LP       10051l
#define     TM_CHXWINLEAVE_SP       10052l
#define     TM_CHXLOW               10060l
#define     TM_CHXHIGH              10061l
#define     TM_CHXINWIN             10062l
#define     TM_CHXOUTWIN            10063l
#define     TM_CHXSPIKE             10064l


#define     TM_CH0POS               10000l
#define     TM_CH0NEG               10010l
#define     TM_CH0OFF               10020l
#define     TM_CH0BOTH              10030l
#define     TM_CH1POS               10100l
#define     TM_CH1NEG               10110l
#define     TM_CH1OFF               10120l
#define     TM_CH1BOTH              10130l
#define     TM_CH2POS               10200l
#define     TM_CH2NEG               10210l
#define     TM_CH2OFF               10220l
#define     TM_CH2BOTH              10230l
#define     TM_CH3POS               10300l
#define     TM_CH3NEG               10310l
#define     TM_CH3OFF               10320l
#define     TM_CH3BOTH              10330l

#define     TM_TTLPOS               20000l
#define     TM_TTLHIGH_LP           20001l
#define     TM_TTLHIGH_SP           20002l
#define     TM_TTLNEG               20010l
#define     TM_TTLLOW_LP            20011l
#define     TM_TTLLOW_SP            20012l
#define     TM_TTL                  20020l
#define     TM_TTLBOTH              20030l
#define     TM_TTLBOTH_LP           20031l
#define     TM_TTLBOTH_SP           20032l
#define     TM_CHANNEL              20040l
#define     TM_TTLHIGH              20050l
#define     TM_TTLLOW               20051l
#define     TM_PATTERN              21000l
#define     TM_PATTERN_LP           21001l
#define     TM_PATTERN_SP           21002l
#define     TM_PATTERNANDEDGE       22000l
#define     TM_PATTERNANDEDGE_LP    22001l
#define     TM_PATTERNANDEDGE_SP    22002l
#define     TM_GATELOW              30000l
#define     TM_GATEHIGH             30010l
#define     TM_GATEPATTERN          30020l
#define     TM_CHOR                 35000l
#define     TM_CHAND                35010l
#define     TM_CHORTTLPOS           35020l
#define     TM_CHORTTLNEG           35021l

#define SPC_PXITRGOUT               40300l
#define     PTO_OFF                  0l
#define     PTO_LINE0                1l
#define     PTO_LINE1                2l
#define     PTO_LINE2                3l
#define     PTO_LINE3                4l
#define     PTO_LINE4                5l
#define     PTO_LINE5                6l
#define     PTO_LINE6                7l
#define     PTO_LINE7                8l
#define     PTO_LINESTAR             9l
#define SPC_PXITRGOUT_AVAILABLE     40301l  // bitmap register

#define SPC_PXISTARTRG_DIVRST_OUT   40302l  // bitmap register
#define SPC_PXISTARTRG_DIVRST_OUT_AVAILABLE   40303l
#define SPC_PXISTARTRG_OUT          40304l  // bitmap register
#define     PSTO_LINESTAR0           0x00000001
#define     PSTO_LINESTAR1           0x00000002
#define     PSTO_LINESTAR2           0x00000004
#define     PSTO_LINESTAR3           0x00000008
#define     PSTO_LINESTAR4           0x00000010
#define     PSTO_LINESTAR5           0x00000020
#define     PSTO_LINESTAR6           0x00000040
#define     PSTO_LINESTAR7           0x00000080
#define     PSTO_LINESTAR8           0x00000100
#define     PSTO_LINESTAR9           0x00000200
#define     PSTO_LINESTAR10          0x00000400
#define     PSTO_LINESTAR11          0x00000800
#define     PSTO_LINESTAR12          0x00001000
#define     PSTO_LINE0               0x00010000
#define     PSTO_LINE1               0x00020000
#define     PSTO_LINE2               0x00040000
#define     PSTO_LINE3               0x00080000
#define     PSTO_LINE4               0x00100000
#define     PSTO_LINE5               0x00200000
#define     PSTO_LINE6               0x00400000
#define     PSTO_LINE7               0x00800000
#define SPC_PXISTARTRG_OUT_AVAILABLE          40305l

#define SPC_PXITRGIN                40310l  // bitmap register
#define     PTI_OFF                  0l
#define     PTI_LINE0                1l
#define     PTI_LINE1                2l
#define     PTI_LINE2                4l
#define     PTI_LINE3                8l
#define     PTI_LINE4                16l
#define     PTI_LINE5                32l
#define     PTI_LINE6                64l
#define     PTI_LINE7                128l
#define     PTI_LINESTAR             256l
#define SPC_PXITRGIN_AVAILABLE      40311l  // bitmap register
#define SPC_PXI_DIVIDER_RESET_IN            40320l
#define SPC_PXI_DIVIDER_RESET_IN_AVAILABLE  40321l


// new registers of M2i driver
#define SPC_TRIG_AVAILORMASK        40400l
#define SPC_TRIG_ORMASK             40410l
#define SPC_TRIG_AVAILANDMASK       40420l
#define SPC_TRIG_ANDMASK            40430l
#define     SPC_TMASK_NONE          0x00000000
#define     SPC_TMASK_SOFTWARE      0x00000001
#define     SPC_TMASK_EXT0          0x00000002
#define     SPC_TMASK_EXT1          0x00000004
#define     SPC_TMASK_EXT2          0x00000008
#define     SPC_TMASK_EXT3          0x00000010
#define     SPC_TMASK_XIO0          0x00000100
#define     SPC_TMASK_XIO1          0x00000200
#define     SPC_TMASK_XIO2          0x00000400
#define     SPC_TMASK_XIO3          0x00000800
#define     SPC_TMASK_XIO4          0x00001000
#define     SPC_TMASK_XIO5          0x00002000
#define     SPC_TMASK_XIO6          0x00004000
#define     SPC_TMASK_XIO7          0x00008000
#define     SPC_TMASK_PXI0          0x00100000
#define     SPC_TMASK_PXI1          0x00200000
#define     SPC_TMASK_PXI2          0x00400000
#define     SPC_TMASK_PXI3          0x00800000
#define     SPC_TMASK_PXI4          0x01000000
#define     SPC_TMASK_PXI5          0x02000000
#define     SPC_TMASK_PXI6          0x04000000
#define     SPC_TMASK_PXI7          0x08000000
#define     SPC_TMASK_PXISTAR       0x10000000
#define     SPC_TMASK_PXIDSTARB     0x20000000

#define SPC_TRIG_CH_AVAILORMASK0    40450l
#define SPC_TRIG_CH_AVAILORMASK1    40451l
#define SPC_TRIG_CH_ORMASK0         40460l
#define SPC_TRIG_CH_ORMASK1         40461l
#define SPC_TRIG_CH_AVAILANDMASK0   40470l
#define SPC_TRIG_CH_AVAILANDMASK1   40471l
#define SPC_TRIG_CH_ANDMASK0        40480l 
#define SPC_TRIG_CH_ANDMASK1        40481l 
#define     SPC_TMASK0_NONE         0x00000000
#define     SPC_TMASK0_CH0          0x00000001
#define     SPC_TMASK0_CH1          0x00000002
#define     SPC_TMASK0_CH2          0x00000004
#define     SPC_TMASK0_CH3          0x00000008
#define     SPC_TMASK0_CH4          0x00000010
#define     SPC_TMASK0_CH5          0x00000020
#define     SPC_TMASK0_CH6          0x00000040
#define     SPC_TMASK0_CH7          0x00000080
#define     SPC_TMASK0_CH8          0x00000100
#define     SPC_TMASK0_CH9          0x00000200
#define     SPC_TMASK0_CH10         0x00000400
#define     SPC_TMASK0_CH11         0x00000800
#define     SPC_TMASK0_CH12         0x00001000
#define     SPC_TMASK0_CH13         0x00002000
#define     SPC_TMASK0_CH14         0x00004000
#define     SPC_TMASK0_CH15         0x00008000
#define     SPC_TMASK0_CH16         0x00010000
#define     SPC_TMASK0_CH17         0x00020000
#define     SPC_TMASK0_CH18         0x00040000
#define     SPC_TMASK0_CH19         0x00080000
#define     SPC_TMASK0_CH20         0x00100000
#define     SPC_TMASK0_CH21         0x00200000
#define     SPC_TMASK0_CH22         0x00400000
#define     SPC_TMASK0_CH23         0x00800000
#define     SPC_TMASK0_CH24         0x01000000
#define     SPC_TMASK0_CH25         0x02000000
#define     SPC_TMASK0_CH26         0x04000000
#define     SPC_TMASK0_CH27         0x08000000
#define     SPC_TMASK0_CH28         0x10000000
#define     SPC_TMASK0_CH29         0x20000000
#define     SPC_TMASK0_CH30         0x40000000
#define     SPC_TMASK0_CH31         0x80000000

#define     SPC_TMASK1_NONE         0x00000000
#define     SPC_TMASK1_CH32         0x00000001
#define     SPC_TMASK1_CH33         0x00000002
#define     SPC_TMASK1_CH34         0x00000004
#define     SPC_TMASK1_CH35         0x00000008
#define     SPC_TMASK1_CH36         0x00000010
#define     SPC_TMASK1_CH37         0x00000020
#define     SPC_TMASK1_CH38         0x00000040
#define     SPC_TMASK1_CH39         0x00000080
#define     SPC_TMASK1_CH40         0x00000100
#define     SPC_TMASK1_CH41         0x00000200
#define     SPC_TMASK1_CH42         0x00000400
#define     SPC_TMASK1_CH43         0x00000800
#define     SPC_TMASK1_CH44         0x00001000
#define     SPC_TMASK1_CH45         0x00002000
#define     SPC_TMASK1_CH46         0x00004000
#define     SPC_TMASK1_CH47         0x00008000
#define     SPC_TMASK1_CH48         0x00010000
#define     SPC_TMASK1_CH49         0x00020000
#define     SPC_TMASK1_CH50         0x00040000
#define     SPC_TMASK1_CH51         0x00080000
#define     SPC_TMASK1_CH52         0x00100000
#define     SPC_TMASK1_CH53         0x00200000
#define     SPC_TMASK1_CH54         0x00400000
#define     SPC_TMASK1_CH55         0x00800000
#define     SPC_TMASK1_CH56         0x01000000
#define     SPC_TMASK1_CH57         0x02000000
#define     SPC_TMASK1_CH58         0x04000000
#define     SPC_TMASK1_CH59         0x08000000
#define     SPC_TMASK1_CH60         0x10000000
#define     SPC_TMASK1_CH61         0x20000000
#define     SPC_TMASK1_CH62         0x40000000
#define     SPC_TMASK1_CH63         0x80000000

#define SPC_TRIG_EXT_AVAILMODES     40500l
#define SPC_TRIG_EXT0_AVAILMODES    40500l
#define SPC_TRIG_EXT1_AVAILMODES    40501l
#define SPC_TRIG_EXT2_AVAILMODES    40502l
#define SPC_TRIG_EXT0_AVAILMODESOR  40503l
#define SPC_TRIG_EXT1_AVAILMODESOR  40504l
#define SPC_TRIG_EXT2_AVAILMODESOR  40505l
#define SPC_TRIG_EXT0_AVAILMODESAND 40506l
#define SPC_TRIG_EXT1_AVAILMODESAND 40507l
#define SPC_TRIG_EXT2_AVAILMODESAND 40508l
#define SPC_TRIG_EXT3_AVAILMODESAND 40509l
#define SPC_TRIG_EXT0_MODE          40510l
#define SPC_TRIG_EXT1_MODE          40511l
#define SPC_TRIG_EXT2_MODE          40512l
#define SPC_TRIG_EXT3_MODE          40513l
#define SPC_TRIG_EXT3_AVAILMODES    40514l
#define SPC_TRIG_EXT3_AVAILMODESOR  40515l

#define SPC_TRIG_EXT0_READFEATURES  40520l
#define SPC_TRIG_EXT1_READFEATURES  40521l
#define SPC_TRIG_EXT2_READFEATURES  40522l
#define SPC_TRIG_EXT3_READFEATURES  40523l
#define     SPCM_TRFEAT_TERM            0x00000001
#define     SPCM_TRFEAT_HIGHIMP         0x00000002
#define     SPCM_TRFEAT_DCCOUPLING      0x00000004
#define     SPCM_TRFEAT_ACCOUPLING      0x00000008
#define     SPCM_TRFEAT_SE              0x00000010
#define     SPCM_TRFEAT_DIFF            0x00000020
#define     SPCM_TRFEAT_LEVELPROG       0x00000100
#define     SPCM_TRFEAT_PROGTHRESHOLD   0x00000200

// legacy constants: not enough contiguous constants possible for X4..X19
#define SPC_LEGACY_X0_READFEATURES  40530l
#define SPC_LEGACY_X1_READFEATURES  40531l
#define SPC_LEGACY_X2_READFEATURES  40532l
#define SPC_LEGACY_X3_READFEATURES  40533l

// legacy constants: not enough contiguous constants possible for X4..X19
#define SPC_LEGACY_X0_TERM          40535l
#define SPC_LEGACY_X1_TERM          40536l
#define SPC_LEGACY_X2_TERM          40537l
#define SPC_LEGACY_X3_TERM          40538l

#define SPC_TRIG_XIO_AVAILMODES     40550l
#define SPC_TRIG_XIO_AVAILMODESOR   40551l
#define SPC_TRIG_XIO_AVAILMODESAND  40552l
#define SPC_TRIG_XIO0_MODE          40560l
#define SPC_TRIG_XIO1_MODE          40561l
#define     SPC_TM_MODEMASK         0x00FFFFFF
#define     SPC_TM_NONE             0x00000000
#define     SPC_TM_POS              0x00000001
#define     SPC_TM_NEG              0x00000002
#define     SPC_TM_BOTH             0x00000004
#define     SPC_TM_HIGH             0x00000008
#define     SPC_TM_LOW              0x00000010
#define     SPC_TM_WINENTER         0x00000020
#define     SPC_TM_WINLEAVE         0x00000040
#define     SPC_TM_INWIN            0x00000080
#define     SPC_TM_OUTSIDEWIN       0x00000100
#define     SPC_TM_SPIKE            0x00000200
#define     SPC_TM_PATTERN          0x00000400
#define     SPC_TM_STEEPPOS         0x00000800
#define     SPC_TM_STEEPNEG         0x00001000
#define     SPC_TM_EXTRAMASK        0xFF000000
#define     SPC_TM_REARM            0x01000000
#define     SPC_TM_PW_SMALLER       0x02000000
#define     SPC_TM_PW_GREATER       0x04000000
#define     SPC_TM_DOUBLEEDGE       0x08000000
#define     SPC_TM_PULSESTRETCH     0x10000000
#define     SPC_TM_HYSTERESIS       0x20000000

#define SPC_TRIG_PATTERN_AVAILMODES 40580l
#define SPC_TRIG_PATTERN_MODE       40590l

#define SPC_TRIG_CH_AVAILMODES      40600l
#define SPC_TRIG_CH_AVAILMODESOR    40601l
#define SPC_TRIG_CH_AVAILMODESAND   40602l
#define SPC_TRIG_CH0_MODE           40610l
#define SPC_TRIG_CH1_MODE           40611l
#define SPC_TRIG_CH2_MODE           40612l
#define SPC_TRIG_CH3_MODE           40613l
#define SPC_TRIG_CH4_MODE           40614l
#define SPC_TRIG_CH5_MODE           40615l
#define SPC_TRIG_CH6_MODE           40616l
#define SPC_TRIG_CH7_MODE           40617l
#define SPC_TRIG_CH8_MODE           40618l
#define SPC_TRIG_CH9_MODE           40619l
#define SPC_TRIG_CH10_MODE          40620l
#define SPC_TRIG_CH11_MODE          40621l
#define SPC_TRIG_CH12_MODE          40622l
#define SPC_TRIG_CH13_MODE          40623l
#define SPC_TRIG_CH14_MODE          40624l
#define SPC_TRIG_CH15_MODE          40625l
#define SPC_TRIG_CH16_MODE          40626l
#define SPC_TRIG_CH17_MODE          40627l
#define SPC_TRIG_CH18_MODE          40628l
#define SPC_TRIG_CH19_MODE          40629l
#define SPC_TRIG_CH20_MODE          40630l
#define SPC_TRIG_CH21_MODE          40631l
#define SPC_TRIG_CH22_MODE          40632l
#define SPC_TRIG_CH23_MODE          40633l
#define SPC_TRIG_CH24_MODE          40634l
#define SPC_TRIG_CH25_MODE          40635l
#define SPC_TRIG_CH26_MODE          40636l
#define SPC_TRIG_CH27_MODE          40637l
#define SPC_TRIG_CH28_MODE          40638l
#define SPC_TRIG_CH29_MODE          40639l
#define SPC_TRIG_CH30_MODE          40640l
#define SPC_TRIG_CH31_MODE          40641l

#define SPC_TRIG_CH32_MODE          40642l
#define SPC_TRIG_CH33_MODE          40643l
#define SPC_TRIG_CH34_MODE          40644l
#define SPC_TRIG_CH35_MODE          40645l
#define SPC_TRIG_CH36_MODE          40646l
#define SPC_TRIG_CH37_MODE          40647l
#define SPC_TRIG_CH38_MODE          40648l
#define SPC_TRIG_CH39_MODE          40649l
#define SPC_TRIG_CH40_MODE          40650l
#define SPC_TRIG_CH41_MODE          40651l
#define SPC_TRIG_CH42_MODE          40652l
#define SPC_TRIG_CH43_MODE          40653l
#define SPC_TRIG_CH44_MODE          40654l
#define SPC_TRIG_CH45_MODE          40655l
#define SPC_TRIG_CH46_MODE          40656l
#define SPC_TRIG_CH47_MODE          40657l
#define SPC_TRIG_CH48_MODE          40658l
#define SPC_TRIG_CH49_MODE          40659l
#define SPC_TRIG_CH50_MODE          40660l
#define SPC_TRIG_CH51_MODE          40661l
#define SPC_TRIG_CH52_MODE          40662l
#define SPC_TRIG_CH53_MODE          40663l
#define SPC_TRIG_CH54_MODE          40664l
#define SPC_TRIG_CH55_MODE          40665l
#define SPC_TRIG_CH56_MODE          40666l
#define SPC_TRIG_CH57_MODE          40667l
#define SPC_TRIG_CH58_MODE          40668l
#define SPC_TRIG_CH59_MODE          40669l
#define SPC_TRIG_CH60_MODE          40670l
#define SPC_TRIG_CH61_MODE          40671l
#define SPC_TRIG_CH62_MODE          40672l
#define SPC_TRIG_CH63_MODE          40673l


#define SPC_TRIG_AVAILDELAY         40800l
#define SPC_TRIG_AVAILDELAY_STEP    40801l
#define SPC_TRIG_DELAY              40810l

#define SPC_TRIG_AVAILHOLDOFF       40802l
#define SPC_TRIG_AVAILHOLDOFF_STEP  40803l
#define SPC_TRIG_HOLDOFF            40811l

#define SPC_SINGLESHOT              41000l
#define SPC_OUTONTRIGGER            41100l
#define SPC_RESTARTCONT             41200l
#define SPC_SINGLERESTART           41300l

#define SPC_TRIGGERLEVEL            42000l
#define SPC_TRIGGERLEVEL0           42000l
#define SPC_TRIGGERLEVEL1           42001l
#define SPC_TRIGGERLEVEL2           42002l
#define SPC_TRIGGERLEVEL3           42003l
#define SPC_TRIGGERLEVEL4           42004l
#define SPC_TRIGGERLEVEL5           42005l
#define SPC_TRIGGERLEVEL6           42006l
#define SPC_TRIGGERLEVEL7           42007l
#define SPC_TRIGGERLEVEL8           42008l
#define SPC_TRIGGERLEVEL9           42009l
#define SPC_TRIGGERLEVEL10          42010l
#define SPC_TRIGGERLEVEL11          42011l
#define SPC_TRIGGERLEVEL12          42012l
#define SPC_TRIGGERLEVEL13          42013l
#define SPC_TRIGGERLEVEL14          42014l
#define SPC_TRIGGERLEVEL15          42015l

#define SPC_AVAILHIGHLEVEL_MIN      41997l
#define SPC_AVAILHIGHLEVEL_MAX      41998l
#define SPC_AVAILHIGHLEVEL_STEP     41999l

#define SPC_HIGHLEVEL0              42000l
#define SPC_HIGHLEVEL1              42001l
#define SPC_HIGHLEVEL2              42002l
#define SPC_HIGHLEVEL3              42003l
#define SPC_HIGHLEVEL4              42004l
#define SPC_HIGHLEVEL5              42005l
#define SPC_HIGHLEVEL6              42006l
#define SPC_HIGHLEVEL7              42007l
#define SPC_HIGHLEVEL8              42008l
#define SPC_HIGHLEVEL9              42009l
#define SPC_HIGHLEVEL10             42010l
#define SPC_HIGHLEVEL11             42011l
#define SPC_HIGHLEVEL12             42012l
#define SPC_HIGHLEVEL13             42013l
#define SPC_HIGHLEVEL14             42014l
#define SPC_HIGHLEVEL15             42015l

#define SPC_AVAILLOWLEVEL_MIN       42097l
#define SPC_AVAILLOWLEVEL_MAX       42098l
#define SPC_AVAILLOWLEVEL_STEP      42099l

#define SPC_LOWLEVEL0               42100l
#define SPC_LOWLEVEL1               42101l
#define SPC_LOWLEVEL2               42102l
#define SPC_LOWLEVEL3               42103l
#define SPC_LOWLEVEL4               42104l
#define SPC_LOWLEVEL5               42105l
#define SPC_LOWLEVEL6               42106l
#define SPC_LOWLEVEL7               42107l
#define SPC_LOWLEVEL8               42108l
#define SPC_LOWLEVEL9               42109l
#define SPC_LOWLEVEL10              42110l
#define SPC_LOWLEVEL11              42111l
#define SPC_LOWLEVEL12              42112l
#define SPC_LOWLEVEL13              42113l
#define SPC_LOWLEVEL14              42114l
#define SPC_LOWLEVEL15              42115l

#define SPC_TRIG_CH0_LEVEL0         42200l
#define SPC_TRIG_CH1_LEVEL0         42201l
#define SPC_TRIG_CH2_LEVEL0         42202l
#define SPC_TRIG_CH3_LEVEL0         42203l
#define SPC_TRIG_CH4_LEVEL0         42204l
#define SPC_TRIG_CH5_LEVEL0         42205l
#define SPC_TRIG_CH6_LEVEL0         42206l
#define SPC_TRIG_CH7_LEVEL0         42207l
#define SPC_TRIG_CH8_LEVEL0         42208l
#define SPC_TRIG_CH9_LEVEL0         42209l
#define SPC_TRIG_CH10_LEVEL0        42210l
#define SPC_TRIG_CH11_LEVEL0        42211l
#define SPC_TRIG_CH12_LEVEL0        42212l
#define SPC_TRIG_CH13_LEVEL0        42213l
#define SPC_TRIG_CH14_LEVEL0        42214l
#define SPC_TRIG_CH15_LEVEL0        42215l

#define SPC_TRIG_CH0_LEVEL1         42300l
#define SPC_TRIG_CH1_LEVEL1         42301l
#define SPC_TRIG_CH2_LEVEL1         42302l
#define SPC_TRIG_CH3_LEVEL1         42303l
#define SPC_TRIG_CH4_LEVEL1         42304l
#define SPC_TRIG_CH5_LEVEL1         42305l
#define SPC_TRIG_CH6_LEVEL1         42306l
#define SPC_TRIG_CH7_LEVEL1         42307l
#define SPC_TRIG_CH8_LEVEL1         42308l
#define SPC_TRIG_CH9_LEVEL1         42309l
#define SPC_TRIG_CH10_LEVEL1        42310l
#define SPC_TRIG_CH11_LEVEL1        42311l
#define SPC_TRIG_CH12_LEVEL1        42312l
#define SPC_TRIG_CH13_LEVEL1        42313l
#define SPC_TRIG_CH14_LEVEL1        42314l
#define SPC_TRIG_CH15_LEVEL1        42315l

#define SPC_TRIG_EXT0_LEVEL0        42320l
#define SPC_TRIG_EXT1_LEVEL0        42321l
#define SPC_TRIG_EXT2_LEVEL0        42322l

#define SPC_TRIG_EXT0_LEVEL1        42330l
#define SPC_TRIG_EXT1_LEVEL1        42331l
#define SPC_TRIG_EXT2_LEVEL1        42332l

#define SPC_TRIG_EXT_AVAIL0_MIN     42340l
#define SPC_TRIG_EXT_AVAIL0_MAX     42341l
#define SPC_TRIG_EXT_AVAIL0_STEP    42342l

#define SPC_TRIG_EXT_AVAIL1_MIN     42345l
#define SPC_TRIG_EXT_AVAIL1_MAX     42346l
#define SPC_TRIG_EXT_AVAIL1_STEP    42347l

// threshold levels (for 77xx)
#define SPC_THRESHOLD0              42400l  // threshold level for channel group 0
#define SPC_THRESHOLD1              42401l  // threshold level for channel group 1
#define SPC_THRESHOLD2              42402l  // threshold level for channel group 2
#define SPC_THRESHOLD3              42403l  // threshold level for channel group 3
#define SPC_CLOCK_THRESHOLD         42410l  // threshold level for clock input
#define SPC_TRIG_THRESHOLD          42411l  // threshold level for trigger input
#define SPC_X0X1_THRESHOLD          42412l  // threshold level for X0/X1 input
#define SPC_STROBE_THRESHOLD        42413l  // threshold level for strobe input

#define SPC_AVAILTHRESHOLD_MIN      42420l
#define SPC_AVAILTHRESHOLD_MAX      42421l
#define SPC_AVAILTHRESHOLD_STEP     42422l

#define SPC_CLOCK_AVAILTHRESHOLD_MIN  42423l
#define SPC_CLOCK_AVAILTHRESHOLD_MAX  42424l
#define SPC_CLOCK_AVAILTHRESHOLD_STEP 42425l

#define SPC_TRIG_AVAILTHRESHOLD_MIN  42426l
#define SPC_TRIG_AVAILTHRESHOLD_MAX  42427l
#define SPC_TRIG_AVAILTHRESHOLD_STEP 42428l

#define SPC_TRIGGERPATTERN          43000l
#define SPC_TRIGGERPATTERN0         43000l
#define SPC_TRIGGERPATTERN1         43001l
#define SPC_TRIGGERMASK             43100l
#define SPC_TRIGGERMASK0            43100l
#define SPC_TRIGGERMASK1            43101l

#define SPC_PULSEWIDTH              44000l
#define SPC_PULSEWIDTH0             44000l
#define SPC_PULSEWIDTH1             44001l

#define SPC_TRIG_CH_AVAILPULSEWIDTH 44100l
#define SPC_TRIG_CH_PULSEWIDTH      44101l
#define SPC_TRIG_CH0_PULSEWIDTH     44101l
#define SPC_TRIG_CH1_PULSEWIDTH     44102l
#define SPC_TRIG_CH2_PULSEWIDTH     44103l
#define SPC_TRIG_CH3_PULSEWIDTH     44104l
#define SPC_TRIG_CH4_PULSEWIDTH     44105l
#define SPC_TRIG_CH5_PULSEWIDTH     44106l
#define SPC_TRIG_CH6_PULSEWIDTH     44107l
#define SPC_TRIG_CH7_PULSEWIDTH     44108l
#define SPC_TRIG_CH8_PULSEWIDTH     44109l
#define SPC_TRIG_CH9_PULSEWIDTH     44110l
#define SPC_TRIG_CH10_PULSEWIDTH    44111l
#define SPC_TRIG_CH11_PULSEWIDTH    44112l
#define SPC_TRIG_CH12_PULSEWIDTH    44113l
#define SPC_TRIG_CH13_PULSEWIDTH    44114l
#define SPC_TRIG_CH14_PULSEWIDTH    44115l
#define SPC_TRIG_CH15_PULSEWIDTH    44116l

#define SPC_TRIG_EXT_AVAILPULSEWIDTH 44200l
#define SPC_TRIG_EXT0_PULSEWIDTH    44210l
#define SPC_TRIG_EXT1_PULSEWIDTH    44211l
#define SPC_TRIG_EXT2_PULSEWIDTH    44212l
#define SPC_TRIG_EXT3_PULSEWIDTH    44213l

// available dividers for MICX
#define SPC_READCLOCKDIVCOUNT       44300l
#define SPC_CLOCKDIV0               44301l
#define SPC_CLOCKDIV1               44302l
#define SPC_CLOCKDIV2               44303l
#define SPC_CLOCKDIV3               44304l
#define SPC_CLOCKDIV4               44305l
#define SPC_CLOCKDIV5               44306l
#define SPC_CLOCKDIV6               44307l
#define SPC_CLOCKDIV7               44308l
#define SPC_CLOCKDIV8               44309l
#define SPC_CLOCKDIV9               44310l
#define SPC_CLOCKDIV10              44311l
#define SPC_CLOCKDIV11              44312l
#define SPC_CLOCKDIV12              44313l
#define SPC_CLOCKDIV13              44314l
#define SPC_CLOCKDIV14              44315l
#define SPC_CLOCKDIV15              44316l
#define SPC_CLOCKDIV16              44317l

#define SPC_READTROFFSET            45000l
#define SPC_TRIGGEREDGE             46000l
#define SPC_TRIGGEREDGE0            46000l
#define SPC_TRIGGEREDGE1            46001l
#define     TE_POS                  10000l
#define     TE_NEG                  10010l
#define     TE_BOTH                 10020l
#define     TE_NONE                 10030l


// ----- Timestamp -----
#define CH_TIMESTAMP                9999l

#define SPC_TIMESTAMP_CMD           47000l
#define     TS_RESET                    0l
#define     TS_MODE_DISABLE             10l
#define     TS_MODE_STARTRESET          11l
#define     TS_MODE_STANDARD            12l
#define     TS_MODE_REFCLOCK            13l
#define     TS_MODE_TEST5555            90l
#define     TS_MODE_TESTAAAA            91l
#define     TS_MODE_ZHTEST              92l

// ----- modes for M2i, M3i, M4i, M4x, M2p hardware (bitmap) -----
#define SPC_TIMESTAMP_AVAILMODES    47001l
#define     SPC_TSMODE_DISABLE      0x00000000
#define     SPC_TS_RESET            0x00000001
#define     SPC_TSMODE_STANDARD     0x00000002
#define     SPC_TSMODE_STARTRESET   0x00000004
#define     SPC_TS_RESET_WAITREFCLK 0x00000008
#define     SPC_TSCNT_INTERNAL      0x00000100
#define     SPC_TSCNT_REFCLOCKPOS   0x00000200
#define     SPC_TSCNT_REFCLOCKNEG   0x00000400
#define     SPC_TSFEAT_NONE         0x00000000
#define     SPC_TSFEAT_STORE1STABA  0x00010000
#define     SPC_TSFEAT_INCRMODE     0x00020000
#define     SPC_TSFEAT_INCRMODE12   0x00040000
#define     SPC_TSFEAT_TRGSRC       0x00080000

#define     SPC_TSXIOACQ_DISABLE    0x00000000
#define     SPC_TSXIOACQ_ENABLE     0x00001000
#define     SPC_TSXIOINC_ENABLE     0x00002000
#define     SPC_TSXIOINC12_ENABLE   0x00004000

#define     SPC_TSMODE_MASK         0x000000FF
#define     SPC_TSCNT_MASK          0x00000F00
#define     SPC_TSFEAT_MASK         0x000F0000

#define     SPC_TRGSRC_MASK_CH0       0x00000001
#define     SPC_TRGSRC_MASK_CH1       0x00000002
#define     SPC_TRGSRC_MASK_CH2       0x00000004
#define     SPC_TRGSRC_MASK_CH3       0x00000008
#define     SPC_TRGSRC_MASK_CH4       0x00000010
#define     SPC_TRGSRC_MASK_CH5       0x00000020
#define     SPC_TRGSRC_MASK_CH6       0x00000040
#define     SPC_TRGSRC_MASK_CH7       0x00000080
#define     SPC_TRGSRC_MASK_EXT0      0x00000100
#define     SPC_TRGSRC_MASK_EXT1      0x00000200
#define     SPC_TRGSRC_MASK_FORCE     0x00000400
// space for digital channels using TSXIOACQ_ENABLE of standard multi-purpose lines
#define     SPC_TRGSRC_MASK_PXI0      0x00010000
#define     SPC_TRGSRC_MASK_PXI1      0x00020000
#define     SPC_TRGSRC_MASK_PXI2      0x00040000
#define     SPC_TRGSRC_MASK_PXI3      0x00080000
#define     SPC_TRGSRC_MASK_PXI4      0x00100000
#define     SPC_TRGSRC_MASK_PXI5      0x00200000
#define     SPC_TRGSRC_MASK_PXI6      0x00400000
#define     SPC_TRGSRC_MASK_PXI7      0x00800000
#define     SPC_TRGSRC_MASK_PXISTAR   0x01000000
#define     SPC_TRGSRC_MASK_PXIDSTARB 0x02000000
#define     SPC_TRGSRC_MASK_X0        0x10000000
#define     SPC_TRGSRC_MASK_X1        0x20000000
#define     SPC_TRGSRC_MASK_X2        0x40000000
#define     SPC_TRGSRC_MASK_X3        0x80000000
// space for more digital channels using TSXIOACQ_ENABLE of additional multi-purpose lines (optional)


#define SPC_TIMESTAMP_STATUS        47010l
#define     TS_FIFO_EMPTY               0l
#define     TS_FIFO_LESSHALF            1l
#define     TS_FIFO_MOREHALF            2l
#define     TS_FIFO_OVERFLOW            3l

#define SPC_TIMESTAMP_COUNT         47020l
#define SPC_TIMESTAMP_STARTTIME     47030l
#define SPC_TIMESTAMP_STARTDATE     47031l
#define SPC_TIMESTAMP_FIFO          47040l
#define SPC_TIMESTAMP_TIMEOUT       47045l

#define SPC_TIMESTAMP_RESETMODE     47050l
#define     TS_RESET_POS               10l
#define     TS_RESET_NEG               20l



// ----- Extra I/O module -----
#define SPC_XIO_DIRECTION           47100l
#define     XD_CH0_INPUT                0l
#define     XD_CH0_OUTPUT               1l
#define     XD_CH1_INPUT                0l
#define     XD_CH1_OUTPUT               2l
#define     XD_CH2_INPUT                0l
#define     XD_CH2_OUTPUT               4l
#define SPC_XIO_DIGITALIO           47110l
#define SPC_XIO_ANALOGOUT0          47120l
#define SPC_XIO_ANALOGOUT1          47121l
#define SPC_XIO_ANALOGOUT2          47122l
#define SPC_XIO_ANALOGOUT3          47123l
#define SPC_XIO_WRITEDACS           47130l



// ----- M3i        multi purpose lines (X0, X1        ) 
// ----- M4i + M4x  multi purpose lines (X0, X1, X2    ) 
// ----- M2p        multi purpose lines (X0, X1, X2, X3) and with installed option also (X4 .. X19)

// legacy constants: not enough contiguous constants possible for X4..X19,
// hence new constants for X-modes (SPCM_X0_MODE.. SPCM_X19_MODE) exist further below
#define SPCM_LEGACY_X0_MODE         47200l
#define SPCM_LEGACY_X1_MODE         47201l
#define SPCM_LEGACY_X2_MODE         47202l
#define SPCM_LEGACY_X3_MODE         47203l
#define SPCM_LEGACY_X0_AVAILMODES   47210l
#define SPCM_LEGACY_X1_AVAILMODES   47211l
#define SPCM_LEGACY_X2_AVAILMODES   47212l
#define SPCM_LEGACY_X3_AVAILMODES   47213l
#define     SPCM_XMODE_DISABLE           0x00000000
#define     SPCM_XMODE_ASYNCIN           0x00000001  // used as asynchronous input
#define     SPCM_XMODE_ASYNCOUT          0x00000002  // used as asynchronous output
#define     SPCM_XMODE_DIGIN             0x00000004  // used as synchronous digital input
#define     SPCM_XMODE_DIGOUT            0x00000008  // used as synchronous digital output
#define     SPCM_XMODE_TRIGIN            0x00000010  // used as trigger input
#define     SPCM_XMODE_TRIGOUT           0x00000020  // used as trigger output
#define     SPCM_XMODE_OVROUT            0x00000040  // used as ADC overrange output
#define     SPCM_XMODE_DIGIN2BIT         0x00000080  // used as synchronous digital input, 2bits per channel
#define     SPCM_XMODE_RUNSTATE          0x00000100  // shows the run state of the card (high = run)
#define     SPCM_XMODE_ARMSTATE          0x00000200  // shows the arm state (high = armed for trigger of one single card)
#define     SPCM_XMODE_DIRECTTRIGOUT     0x00000400  // used as direct trigger output (safe mode) 
#define     SPCM_XMODE_DIRECTTRIGOUT_LR  0x00000800  // used as direct trigger output (low re-arm)
#define     SPCM_XMODE_REFCLKOUT         0x00001000  // outputs internal or fed in external refclock
#define     SPCM_XMODE_CONTOUTMARK       0x00002000  // outputs a half posttrigger long HIGH pulse on replay
#define     SPCM_XMODE_SYSCLKOUT         0x00004000  // outputs internal system clock
#define     SPCM_XMODE_CLKOUT            0x00008000  // clock output
#define     SPCM_XMODE_SYNCARMSTATE      0x00010000  // shows the arm state (high = armed for trigger when all cards connected to a Star-Hub are armed)
#define     SPCM_XMODE_OPTDIGIN2BIT      0x00020000  // used as synchronous digital input from digitaloption, 2bits per channel
#define     SPCM_XMODE_OPTDIGIN4BIT      0x00040000  // used as synchronous digital input from digitaloption, 4bits per channel
#define     SPCM_XMODE_MODEMASK          0x000FFFFF

// additional constants to be combined together with SPCM_XMODE_DIGOUT to select analog channel containing digital data
#define     SPCM_XMODE_DIGOUTSRC_CH0     0x01000000  // Select Ch0 as source 
#define     SPCM_XMODE_DIGOUTSRC_CH1     0x02000000  // Select Ch1 as source
#define     SPCM_XMODE_DIGOUTSRC_CH2     0x04000000  // Select Ch2 as source
#define     SPCM_XMODE_DIGOUTSRC_CH3     0x08000000  // Select Ch3 as source
#define     SPCM_XMODE_DIGOUTSRC_CH4     0x10000000  // Select Ch4 as source
#define     SPCM_XMODE_DIGOUTSRC_CH5     0x20000000  // Select Ch5 as source
#define     SPCM_XMODE_DIGOUTSRC_CH6     0x40000000  // Select Ch6 as source
#define     SPCM_XMODE_DIGOUTSRC_CH7     0x80000000  // Select Ch7 as source
#define     SPCM_XMODE_DIGOUTSRC_CHMASK  0xFF000000

// additional constants to be combined together with SPCM_XMODE_DIGOUT to select digital signal source
#define     SPCM_XMODE_DIGOUTSRC_BIT15              0x00100000  // Use Bit15 (MSB    ) of selected channel: channel resolution will be reduced to 15 bit
#define     SPCM_XMODE_DIGOUTSRC_BIT14              0x00200000  // Use Bit14 (MSB - 1) of selected channel: channel resolution will be reduced to 14 bit
#define     SPCM_XMODE_DIGOUTSRC_BIT13              0x00400000  // Use Bit13 (MSB - 2) of selected channel: channel resolution will be reduced to 13 bit
#define     SPCM_XMODE_DIGOUTSRC_BIT12              0x00800000  // Use Bit12 (MSB - 3) of selected channel: channel resolution will be reduced to 12 bit
#define     SPCM_XMODE_DIGOUTSRC_BITMASK            0x00F00000
// special combinations for M2p.65xx cards with options SPCM_FEAT_DIG16_SMB or SPCM_FEAT_DIG16_FX2
#define     SPCM_XMODE_DIGOUTSRC_BIT15_downto_0     0x00F00000  // use all   16 bits of selected channel on  (X19..X4)              : channel will only contain digital data
#define     SPCM_XMODE_DIGOUTSRC_BIT15_downto_8     0x00700000  // use upper  8 bits of selected channel for (X19..X12) or (X11..X4): channel resolution will be reduced to 8 bit

#define SPCM_XX_ASYNCIO             47220l           // asynchronous in/out register

#define SPC_DIGMODE0 47250l
#define SPC_DIGMODE1 47251l
#define SPC_DIGMODE2 47252l
#define SPC_DIGMODE3 47253l
#define SPC_DIGMODE4 47254l
#define SPC_DIGMODE5 47255l
#define SPC_DIGMODE6 47256l
#define SPC_DIGMODE7 47257l
#define     SPCM_DIGMODE_OFF 0x00000000

#define     SPCM_DIGMODE_X1  0x294A5000 // (M2P_DIGMODE_X1 << (32 - 5)) | (M2P_DIGMODE_X1 << (32 - 10))  ... etc
#define     SPCM_DIGMODE_X2  0x318C6000 // (M2P_DIGMODE_X2 << (32 - 5)) | (M2P_DIGMODE_X2 << (32 - 10))  ... etc
#define     SPCM_DIGMODE_X3  0x39CE7000 // (M2P_DIGMODE_X3 << (32 - 5)) | (M2P_DIGMODE_X3 << (32 - 10))  ... etc
#define     SPCM_DIGMODE_X4  0x84210001
#define     SPCM_DIGMODE_X5  0x8c631002
#define     SPCM_DIGMODE_X6  0x94a52004
#define     SPCM_DIGMODE_X7  0x9ce73008
#define     SPCM_DIGMODE_X8  0xa5294010
#define     SPCM_DIGMODE_X9  0xad6b5020
#define     SPCM_DIGMODE_X10 0xb5ad6040
#define     SPCM_DIGMODE_X11 0xbdef7080
#define     SPCM_DIGMODE_X12 0xc6318100
#define     SPCM_DIGMODE_X13 0xce739200
#define     SPCM_DIGMODE_X14 0xd6b5a400
#define     SPCM_DIGMODE_X15 0xdef7b800
#define     SPCM_DIGMODE_X16 0xe739c000
#define     SPCM_DIGMODE_X17 0xef7bd000
#define     SPCM_DIGMODE_X18 0xf7bde000
#define     SPCM_DIGMODE_X19 0xfffff000

#define     DIGMODEMASK_BIT15 0xF8000000
#define     DIGMODEMASK_BIT14 0x07C00000
#define     DIGMODEMASK_BIT13 0x003E0000
#define     DIGMODEMASK_BIT12 0x0001F000
#define     DIGMODEMASK_BIT11 0x00000800 // one bit only for bit 11 downto 0
#define     DIGMODEMASK_BIT10 0x00000400
#define     DIGMODEMASK_BIT9  0x00000200
#define     DIGMODEMASK_BIT8  0x00000100
#define     DIGMODEMASK_BIT7  0x00000080
#define     DIGMODEMASK_BIT6  0x00000040
#define     DIGMODEMASK_BIT5  0x00000020
#define     DIGMODEMASK_BIT4  0x00000010
#define     DIGMODEMASK_BIT3  0x00000008
#define     DIGMODEMASK_BIT2  0x00000004
#define     DIGMODEMASK_BIT1  0x00000002
#define     DIGMODEMASK_BIT0  0x00000001

// provided for convenience
#define SPCM_DIGMODE_CHREPLACE     0xFFBBCFFF
//#define SPCM_DIGMODE_CHREPLACE    (  (DIGMODEMASK_BIT15 & SPCM_DIGMODE_X19)
//                                   | (DIGMODEMASK_BIT14 & SPCM_DIGMODE_X18)
//                                   | (DIGMODEMASK_BIT13 & SPCM_DIGMODE_X17)
//                                   | (DIGMODEMASK_BIT12 & SPCM_DIGMODE_X16)
//                                   | (DIGMODEMASK_BIT11 & SPCM_DIGMODE_X15)
//                                   | (DIGMODEMASK_BIT10 & SPCM_DIGMODE_X14)
//                                   | (DIGMODEMASK_BIT9  & SPCM_DIGMODE_X13)
//                                   | (DIGMODEMASK_BIT8  & SPCM_DIGMODE_X12)
//                                   | (DIGMODEMASK_BIT7  & SPCM_DIGMODE_X11)
//                                   | (DIGMODEMASK_BIT6  & SPCM_DIGMODE_X10)
//                                   | (DIGMODEMASK_BIT5  & SPCM_DIGMODE_X9 )
//                                   | (DIGMODEMASK_BIT4  & SPCM_DIGMODE_X8 )
//                                   | (DIGMODEMASK_BIT3  & SPCM_DIGMODE_X7 )
//                                   | (DIGMODEMASK_BIT2  & SPCM_DIGMODE_X6 )
//                                   | (DIGMODEMASK_BIT1  & SPCM_DIGMODE_X5 )
//                                   | (DIGMODEMASK_BIT0  & SPCM_DIGMODE_X4 ) )
//


// ----- M4x PXI Trigger lines -----
#define SPC_PXITRG0_MODE           47300l
#define SPC_PXITRG1_MODE           47301l
#define SPC_PXITRG2_MODE           47302l
#define SPC_PXITRG3_MODE           47303l
#define SPC_PXITRG4_MODE           47304l
#define SPC_PXITRG5_MODE           47305l
#define SPC_PXITRG6_MODE           47306l
#define SPC_PXITRG7_MODE           47307l
#define SPC_PXISTAR_MODE           47308l
#define SPC_PXIDSTARC_MODE         47309l
#define SPC_PXITRG0_AVAILMODES     47310l
#define SPC_PXITRG1_AVAILMODES     47311l
#define SPC_PXITRG2_AVAILMODES     47312l
#define SPC_PXITRG3_AVAILMODES     47313l
#define SPC_PXITRG4_AVAILMODES     47314l
#define SPC_PXITRG5_AVAILMODES     47315l
#define SPC_PXITRG6_AVAILMODES     47316l
#define SPC_PXITRG7_AVAILMODES     47317l
#define SPC_PXISTAR_AVAILMODES     47318l
#define SPC_PXIDSTARC_AVAILMODES   47319l
#define SPC_PXITRG_ASYNCIO         47320l          // asynchronous in/out register
#define     SPCM_PXITRGMODE_DISABLE     0x00000000
#define     SPCM_PXITRGMODE_IN          0x00000001  // used as input
#define     SPCM_PXITRGMODE_ASYNCOUT    0x00000002  // used as asynchronous output
#define     SPCM_PXITRGMODE_RUNSTATE    0x00000004  // shows the run state of the card (high = run)
#define     SPCM_PXITRGMODE_ARMSTATE    0x00000008  // shows the arm state (high = armed for trigger)
#define     SPCM_PXITRGMODE_TRIGOUT     0x00000010  // used as trigger output
#define     SPCM_PXITRGMODE_REFCLKOUT   0x00000020  // outputs PXI refclock (10 MHz)
#define     SPCM_PXITRGMODE_CONTOUTMARK 0x00000040  // outputs a half posttrigger long HIGH pulse on replay


// ----- Star-Hub -----
// 48000 not usable

#define SPC_STARHUB_STATUS          48010l

#define SPC_STARHUB_ROUTE0          48100l  // Routing Information for Test
#define SPC_STARHUB_ROUTE99         48199l  // ...


// Spcm driver (M2i, M3i, M4i, M4x, M2p) sync setup registers
#define SPC_SYNC_READ_SYNCCOUNT     48990l  // number of sync'd cards
#define SPC_SYNC_READ_NUMCONNECTORS 48991l  // number of connectors on starhub

#define SPC_SYNC_READ_CARDIDX0      49000l  // read index of card at location 0 of sync
#define SPC_SYNC_READ_CARDIDX1      49001l  // ...
#define SPC_SYNC_READ_CARDIDX2      49002l  // ...
#define SPC_SYNC_READ_CARDIDX3      49003l  // ...
#define SPC_SYNC_READ_CARDIDX4      49004l  // ...
#define SPC_SYNC_READ_CARDIDX5      49005l  // ...
#define SPC_SYNC_READ_CARDIDX6      49006l  // ...
#define SPC_SYNC_READ_CARDIDX7      49007l  // ...
#define SPC_SYNC_READ_CARDIDX8      49008l  // ...
#define SPC_SYNC_READ_CARDIDX9      49009l  // ...
#define SPC_SYNC_READ_CARDIDX10     49010l  // ...
#define SPC_SYNC_READ_CARDIDX11     49011l  // ...
#define SPC_SYNC_READ_CARDIDX12     49012l  // ...
#define SPC_SYNC_READ_CARDIDX13     49013l  // ...
#define SPC_SYNC_READ_CARDIDX14     49014l  // ...
#define SPC_SYNC_READ_CARDIDX15     49015l  // ...

#define SPC_SYNC_READ_CABLECON0     49100l  // read cable connection of card at location 0 of sync
#define SPC_SYNC_READ_CABLECON1     49101l  // ...
#define SPC_SYNC_READ_CABLECON2     49102l  // ...
#define SPC_SYNC_READ_CABLECON3     49103l  // ...
#define SPC_SYNC_READ_CABLECON4     49104l  // ...
#define SPC_SYNC_READ_CABLECON5     49105l  // ...
#define SPC_SYNC_READ_CABLECON6     49106l  // ...
#define SPC_SYNC_READ_CABLECON7     49107l  // ...
#define SPC_SYNC_READ_CABLECON8     49108l  // ...
#define SPC_SYNC_READ_CABLECON9     49109l  // ...
#define SPC_SYNC_READ_CABLECON10    49110l  // ...
#define SPC_SYNC_READ_CABLECON11    49111l  // ...
#define SPC_SYNC_READ_CABLECON12    49112l  // ...
#define SPC_SYNC_READ_CABLECON13    49113l  // ...
#define SPC_SYNC_READ_CABLECON14    49114l  // ...
#define SPC_SYNC_READ_CABLECON15    49115l  // ...

#define SPC_SYNC_ENABLEMASK         49200l  // synchronisation enable (mask)
#define SPC_SYNC_NOTRIGSYNCMASK     49210l  // trigger disabled for sync (mask)
#define SPC_SYNC_CLKMASK            49220l  // clock master (mask)
#define SPC_SYNC_MODE               49230l  // synchronization mode
#define SPC_AVAILSYNC_MODES         49231l  // available synchronization modes
#define     SPC_SYNC_STANDARD         0x00000001  // starhub uses its own clock and trigger sources
#define     SPC_SYNC_SYSTEMCLOCK      0x00000002  // starhub uses own trigger sources and takes clock from system starhub
#define     SPC_SYNC_SYSTEMCLOCKTRIG  0x00000004  // starhub takes clock and trigger from system starhub (trigger sampled on rising  clock edge)
#define     SPC_SYNC_SYSTEMCLOCKTRIGN 0x00000008  // starhub takes clock and trigger from system starhub (trigger sampled on falling clock edge)
#define SPC_SYNC_SYSTEM_TRIGADJUST  49240l  // Delay value for adjusting trigger position using system starhub


// ----- Gain and Offset Adjust DAC's -----
#define SPC_ADJ_START               50000l

#define SPC_ADJ_LOAD                50000l
#define SPC_ADJ_SAVE                50010l
#define     ADJ_DEFAULT                 0l
#define     ADJ_USER0                   1l
#define     ADJ_USER1                   2l
#define     ADJ_USER2                   3l
#define     ADJ_USER3                   4l
#define     ADJ_USER4                   5l
#define     ADJ_USER5                   6l
#define     ADJ_USER6                   7l
#define     ADJ_USER7                   8l

#define SPC_ADJ_AUTOADJ             50020l
#define     ADJ_ALL                     0l
#define     ADJ_CURRENT                 1l
#define     ADJ_EXTERNAL                2l
#define     ADJ_1MOHM                   3l

#define     ADJ_CURRENT_CLOCK           4l
#define     ADJ_CURRENT_IR              8l
#define     ADJ_OFFSET_ONLY            16l
#define     ADJ_SPECIAL_CLOCK          32l

#define SPC_ADJ_SOURCE_CALLBACK     50021l
#define SPC_ADJ_PROGRESS_CALLBACK   50022l

#define SPC_ADJ_SET                 50030l
#define SPC_ADJ_FAILMASK            50040l

#define SPC_ADJ_CALIBSOURCE            50050l
#define        ADJ_CALSRC_GAIN             1l
#define        ADJ_CALSRC_OFF              0l
#define        ADJ_CALSRC_GND             -1l
#define        ADJ_CALSRC_GNDOFFS         -2l
#define        ADJ_CALSRC_AC              10l

#define SPC_ADJ_CALIBVALUE0            50060l
#define SPC_ADJ_CALIBVALUE1            50061l
#define SPC_ADJ_CALIBVALUE2            50062l
#define SPC_ADJ_CALIBVALUE3            50063l
#define SPC_ADJ_CALIBVALUE4            50064l
#define SPC_ADJ_CALIBVALUE5            50065l
#define SPC_ADJ_CALIBVALUE6            50066l
#define SPC_ADJ_CALIBVALUE7            50067l

#define SPC_ADJ_OFFSET_CH0          50900l
#define SPC_ADJ_OFFSET_CH1          50901l
#define SPC_ADJ_OFFSET_CH2          50902l
#define SPC_ADJ_OFFSET_CH3          50903l
#define SPC_ADJ_OFFSET_CH4          50904l
#define SPC_ADJ_OFFSET_CH5          50905l
#define SPC_ADJ_OFFSET_CH6          50906l
#define SPC_ADJ_OFFSET_CH7          50907l
#define SPC_ADJ_OFFSET_CH8          50908l
#define SPC_ADJ_OFFSET_CH9          50909l
#define SPC_ADJ_OFFSET_CH10         50910l
#define SPC_ADJ_OFFSET_CH11         50911l
#define SPC_ADJ_OFFSET_CH12         50912l
#define SPC_ADJ_OFFSET_CH13         50913l
#define SPC_ADJ_OFFSET_CH14         50914l
#define SPC_ADJ_OFFSET_CH15         50915l

#define SPC_ADJ_GAIN_CH0            50916l
#define SPC_ADJ_GAIN_CH1            50917l
#define SPC_ADJ_GAIN_CH2            50918l
#define SPC_ADJ_GAIN_CH3            50919l
#define SPC_ADJ_GAIN_CH4            50920l
#define SPC_ADJ_GAIN_CH5            50921l
#define SPC_ADJ_GAIN_CH6            50922l
#define SPC_ADJ_GAIN_CH7            50923l
#define SPC_ADJ_GAIN_CH8            50924l
#define SPC_ADJ_GAIN_CH9            50925l
#define SPC_ADJ_GAIN_CH10           50926l
#define SPC_ADJ_GAIN_CH11           50927l
#define SPC_ADJ_GAIN_CH12           50928l
#define SPC_ADJ_GAIN_CH13           50929l
#define SPC_ADJ_GAIN_CH14           50930l
#define SPC_ADJ_GAIN_CH15           50931l

#define SPC_ADJ_OFFSET0             51000l
#define SPC_ADJ_OFFSET999           51999l

#define SPC_ADJ_GAIN0               52000l
#define SPC_ADJ_GAIN999             52999l

#define SPC_ADJ_CORRECT0            53000l
#define SPC_ADJ_OFFS_CORRECT0       53000l
#define SPC_ADJ_CORRECT999          53999l
#define SPC_ADJ_OFFS_CORRECT999     53999l

#define SPC_ADJ_XIOOFFS0            54000l
#define SPC_ADJ_XIOOFFS1            54001l
#define SPC_ADJ_XIOOFFS2            54002l
#define SPC_ADJ_XIOOFFS3            54003l

#define SPC_ADJ_XIOGAIN0            54010l
#define SPC_ADJ_XIOGAIN1            54011l
#define SPC_ADJ_XIOGAIN2            54012l
#define SPC_ADJ_XIOGAIN3            54013l

#define SPC_ADJ_GAIN_CORRECT0       55000l
#define SPC_ADJ_GAIN_CORRECT999     55999l

#define SPC_ADJ_OFFSCALIBCORRECT0   56000l
#define SPC_ADJ_OFFSCALIBCORRECT999 56999l

#define SPC_ADJ_GAINCALIBCORRECT0   57000l
#define SPC_ADJ_GAINCALIBCORRECT999 57999l

#define SPC_ADJ_ANALOGTRIGGER0      58000l
#define SPC_ADJ_ANALOGTRIGGER99     58099l

#define SPC_ADJ_CALIBSAMPLERATE0    58100l
#define SPC_ADJ_CALIBSAMPLERATE99   58199l

#define SPC_ADJ_CALIBSAMPLERATE_GAIN0    58200l
#define SPC_ADJ_CALIBSAMPLERATE_GAIN99   58299l

#define SPC_ADJ_REFCLOCK            58300l
#define SPC_ADJ_STARHUB_REFCLOCK    58301l

#define SPC_ADJ_END                 59999l



// ----- FIFO Control -----
#define SPC_FIFO_BUFFERS            60000l          // number of FIFO buffers
#define SPC_FIFO_BUFLEN             60010l          // len of each FIFO buffer
#define SPC_FIFO_BUFCOUNT           60020l          // number of FIFO buffers tranfered until now
#define SPC_FIFO_BUFMAXCNT          60030l          // number of FIFO buffers to be transfered (0=continuous)
#define SPC_FIFO_BUFADRCNT          60040l          // number of FIFO buffers allowed
#define SPC_FIFO_BUFREADY           60050l          // fifo buffer ready register (same as SPC_COMMAND + SPC_FIFO_BUFREADY0...)
#define SPC_FIFO_BUFFILLCNT         60060l          // number of currently filled buffers
#define SPC_FIFO_BUFADR0            60100l          // adress of FIFO buffer no. 0
#define SPC_FIFO_BUFADR1            60101l          // ...
#define SPC_FIFO_BUFADR2            60102l          // ...
#define SPC_FIFO_BUFADR3            60103l          // ...
#define SPC_FIFO_BUFADR4            60104l          // ...
#define SPC_FIFO_BUFADR5            60105l          // ...
#define SPC_FIFO_BUFADR6            60106l          // ...
#define SPC_FIFO_BUFADR7            60107l          // ...
#define SPC_FIFO_BUFADR8            60108l          // ...
#define SPC_FIFO_BUFADR9            60109l          // ...
#define SPC_FIFO_BUFADR10           60110l          // ...
#define SPC_FIFO_BUFADR11           60111l          // ...
#define SPC_FIFO_BUFADR12           60112l          // ...
#define SPC_FIFO_BUFADR13           60113l          // ...
#define SPC_FIFO_BUFADR14           60114l          // ...
#define SPC_FIFO_BUFADR15           60115l          // ...
#define SPC_FIFO_BUFADR255          60355l          // last



// ----- Filter -----
#define SPC_FILTER                  100000l
#define SPC_READNUMFILTERS          100001l         // number of programable filters
#define SPC_FILTERFREQUENCY0        100002l         // frequency of filter 0 (bypass)
#define SPC_FILTERFREQUENCY1        100003l         // frequency of filter 1
#define SPC_FILTERFREQUENCY2        100004l         // frequency of filter 2
#define SPC_FILTERFREQUENCY3        100005l         // frequency of filter 3
#define SPC_DIGITALBWFILTER         100100l         // enable/disable digital bandwith filter


// ----- Pattern -----
#define SPC_PATTERNENABLE           110000l
#define SPC_READDIGITAL             110100l

#define SPC_DIGITALMODE0            110200l
#define SPC_DIGITALMODE1            110201l
#define SPC_DIGITALMODE2            110202l
#define SPC_DIGITALMODE3            110203l
#define SPC_DIGITALMODE4            110204l
#define SPC_DIGITALMODE5            110205l
#define SPC_DIGITALMODE6            110206l
#define SPC_DIGITALMODE7            110207l
#define     SPC_DIGITALMODE_OFF         0l
#define     SPC_DIGITALMODE_2BIT        1l
#define     SPC_DIGITALMODE_4BIT        2l
#define     SPC_DIGITALMODE_CHREPLACE   3l


// ----- Miscellanous -----
#define SPC_MISCDAC0                200000l
#define SPC_MISCDAC1                200010l
#define SPC_FACTORYMODE             200020l
#define SPC_DIRECTDAC               200030l
#define SPC_NOTRIGSYNC              200040l
#define SPC_DSPDIRECT               200100l
#define SPC_DMAPHYSICALADR          200110l
#define SPC_MICXCOMP_CLOSEBOARD     200119l
#define SPC_MICXCOMPATIBILITYMODE   200120l
#define SPC_TEST_FIFOSPEED          200121l
#define SPC_RELOADDEMO              200122l
#define SPC_OVERSAMPLINGFACTOR      200123l
#define SPC_ISMAPPEDCARD            200124l
#define     SPCM_NOT_MAPPED             0l
#define     SPCM_LOCAL_MAPPED           1l
#define     SPCM_REMOTE_MAPPED          2l
#define SPC_GETTHREADHANDLE         200130l
#define SPC_GETKERNELHANDLE         200131l
#define SPC_XYZMODE                 200200l
#define SPC_INVERTDATA              200300l
#define SPC_GATEMARKENABLE          200400l
#define SPC_GATE_LEN_ALIGNMENT      200401l
#define SPC_CONTOUTMARK             200450l
#define SPC_EXPANDINT32             200500l
#define SPC_NOPRETRIGGER            200600l
#define SPC_RELAISWAITTIME          200700l
#define SPC_DACWAITTIME             200710l
#define SPC_DELAY_US                200720l
#define SPC_ILAMODE                 200800l
#define SPC_NMDGMODE                200810l
#define SPC_CKADHALF_OUTPUT         200820l
#define SPC_LONGTRIG_OUTPUT         200830l
#define SPC_STOREMODAENDOFSEGMENT   200840l
#define SPC_COUNTERMODE             200850l
#define     SPC_CNTMOD_MASK             0x0000000F
#define     SPC_CNTMOD_PARALLELDATA     0x00000000
#define     SPC_CNTMOD_8BITCNT          0x00000001
#define     SPC_CNTMOD_2x8BITCNT        0x00000002
#define     SPC_CNTMOD_16BITCNT         0x00000003
#define     SPC_CNT0_MASK               0x000000F0
#define     SPC_CNT0_CNTONPOSEDGE       0x00000000
#define     SPC_CNT0_CNTONNEGEDGE       0x00000010
#define     SPC_CNT0_RESETHIGHLVL       0x00000000
#define     SPC_CNT0_RESETLOWLVL        0x00000020
#define     SPC_CNT0_STOPATMAX          0x00000000
#define     SPC_CNT0_ROLLOVER           0x00000040
#define     SPC_CNT1_MASK               0x00000F00
#define     SPC_CNT1_CNTONPOSEDGE       0x00000000
#define     SPC_CNT1_CNTONNEGEDGE       0x00000100
#define     SPC_CNT1_RESETHIGHLVL       0x00000000
#define     SPC_CNT1_RESETLOWLVL        0x00000200
#define     SPC_CNT1_STOPATMAX          0x00000000
#define     SPC_CNT1_ROLLOVER           0x00000400
#define     SPC_CNTCMD_MASK             0x0000F000
#define     SPC_CNTCMD_RESETCNT0        0x00001000
#define     SPC_CNTCMD_RESETCNT1        0x00002000
#define SPC_ENHANCEDSTATUS          200900l
#define     SPC_ENHSTAT_OVERRANGE0      0x00000001
#define     SPC_ENHSTAT_OVERRANGE1      0x00000002
#define     SPC_ENHSTAT_OVERRANGE2      0x00000004
#define     SPC_ENHSTAT_OVERRANGE3      0x00000008
#define     SPC_ENHSTAT_OVERRANGE4      0x00000010
#define     SPC_ENHSTAT_OVERRANGE5      0x00000020
#define     SPC_ENHSTAT_OVERRANGE6      0x00000040
#define     SPC_ENHSTAT_OVERRANGE7      0x00000080
#define     SPC_ENHSTAT_COMPARATOR0     0x40000000
#define     SPC_ENHSTAT_COMPARATOR1     0x80000000
#define     SPC_ENHSTAT_COMPARATOR2     0x20000000
#define     SPC_ENHSTAT_TRGCOMPARATOR   0x40000000
#define     SPC_ENHSTAT_CLKCOMPARATOR   0x80000000
#define SPC_TRIGGERCOUNTER          200905l
#define SPC_FILLSIZEPROMILLE        200910l
#define SPC_OVERRANGEBIT            201000l
#define SPC_2CH8BITMODE             201100l
#define SPC_12BITMODE               201200l
#define SPC_HOLDLASTSAMPLE          201300l

#define SPC_DATACONVERSION          201400l
#define SPC_AVAILDATACONVERSION     201401l
#define     SPCM_DC_NONE            0x00000000
#define     SPCM_DC_12BIT_TO_14BIT  0x00000001
#define     SPCM_DC_16BIT_TO_14BIT  0x00000002
#define     SPCM_DC_12BIT_TO_16BIT  0x00000004
#define     SPCM_DC_14BIT_TO_16BIT  0x00000008
#define     SPCM_DC_15BIT_TO_16BIT  0x00000010
#define     SPCM_DC_13BIT_TO_16BIT  0x00000020
#define     SPCM_DC_14BIT_TO_8BIT   0x00000100
#define     SPCM_DC_16BIT_TO_8BIT   0x00000200
#define     SPCM_DC_16BIT_TO_12BIT  0x00000400
#define     SPCM_DC_TO_OFFSETBINARY 0x00000800

#define SPC_CARDIDENTIFICATION      201500l

#define SPC_HANDSHAKE               201600l

#define SPC_CKSYNC0                 202000l
#define SPC_CKSYNC1                 202001l
#define SPC_DISABLEMOD0             203000l
#define SPC_DISABLEMOD1             203010l
#define SPC_ENABLEOVERRANGECHECK    204000l
#define SPC_OVERRANGESTATUS         204010l
#define SPC_BITMODE                 205000l

#define SPC_READBACK                206000l
#define SPC_AVAILSTOPLEVEL          206009l
#define SPC_STOPLEVEL1              206010l
#define SPC_STOPLEVEL0              206020l
#define SPC_CH0_STOPLEVEL           206020l
#define SPC_CH1_STOPLEVEL           206021l
#define SPC_CH2_STOPLEVEL           206022l
#define SPC_CH3_STOPLEVEL           206023l
#define SPC_CH4_STOPLEVEL           206024l
#define SPC_CH5_STOPLEVEL           206025l
#define SPC_CH6_STOPLEVEL           206026l
#define SPC_CH7_STOPLEVEL           206027l
#define     SPCM_STOPLVL_TRISTATE   0x00000001
#define     SPCM_STOPLVL_LOW        0x00000002
#define     SPCM_STOPLVL_HIGH       0x00000004
#define     SPCM_STOPLVL_HOLDLAST   0x00000008
#define     SPCM_STOPLVL_ZERO       0x00000010
#define     SPCM_STOPLVL_CUSTOM     0x00000020

#define SPC_DIFFMODE                206030l
#define SPC_DACADJUST               206040l

#define SPC_CH0_CUSTOM_STOP         206050l
#define SPC_CH1_CUSTOM_STOP         206051l
#define SPC_CH2_CUSTOM_STOP         206052l
#define SPC_CH3_CUSTOM_STOP         206053l
#define SPC_CH4_CUSTOM_STOP         206054l
#define SPC_CH5_CUSTOM_STOP         206055l
#define SPC_CH6_CUSTOM_STOP         206056l
#define SPC_CH7_CUSTOM_STOP         206057l

#define SPC_AMP_MODE                207000l

#define SPCM_FW_CTRL                210000l
#define SPCM_FW_CTRL_GOLDEN         210001l
#define SPCM_FW_CTRL_ACTIVE         210002l
#define SPCM_FW_CLOCK               210010l
#define SPCM_FW_CONFIG              210020l
#define SPCM_FW_MODULEA             210030l
#define SPCM_FW_MODULEB             210031l
#define SPCM_FW_MODULEA_ACTIVE      210032l
#define SPCM_FW_MODULEB_ACTIVE      210033l
#define SPCM_FW_MODEXTRA            210050l
#define SPCM_FW_MODEXTRA_ACTIVE     210052l
#define SPCM_FW_POWER               210060l
#define SPCM_FW_POWER_ACTIVE        210062l

#define SPC_MULTI                   220000l
#define SPC_DOUBLEMEM               220100l
#define SPC_MULTIMEMVALID           220200l
#define SPC_BANK                    220300l
#define SPC_GATE                    220400l
#define SPC_RELOAD                  230000l
#define SPC_USEROUT                 230010l
#define SPC_WRITEUSER0              230100l
#define SPC_WRITEUSER1              230110l
#define SPC_READUSER0               230200l
#define SPC_READUSER1               230210l
#define SPC_MUX                     240000l
#define SPC_ADJADC                  241000l
#define SPC_ADJOFFS0                242000l
#define SPC_ADJOFFS1                243000l
#define SPC_ADJGAIN0                244000l
#define SPC_ADJGAIN1                245000l
#define SPC_READEPROM               250000l
#define SPC_WRITEEPROM              250010l
#define SPC_DIRECTIO                260000l
#define SPC_DIRECT_MODA             260010l
#define SPC_DIRECT_MODB             260020l
#define SPC_DIRECT_EXT0             260030l
#define SPC_DIRECT_EXT1             260031l
#define SPC_DIRECT_EXT2             260032l
#define SPC_DIRECT_EXT3             260033l
#define SPC_DIRECT_EXT4             260034l
#define SPC_DIRECT_EXT5             260035l
#define SPC_DIRECT_EXT6             260036l
#define SPC_DIRECT_EXT7             260037l
#define SPC_MEMTEST                 270000l
#define SPC_NODMA                   275000l
#define SPC_NOCOUNTER               275010l
#define SPC_NOSCATTERGATHER         275020l
#define SPC_USER_RELAIS_OVERWRITE   275030l
#define     SPCM_URO_ENABLE             0x80000000
#define     SPCM_URO_INVERT_10TO1REL    0x00000001
#define SPC_RUNINTENABLE            290000l
#define SPC_XFERBUFSIZE             295000l
#define SPC_CHLX                    295010l
#define SPC_SPECIALCLOCK            295100l
#define SPC_PLL0_ICP                295105l
#define     SPCM_ICP0            0x00000000
// ...
#define     SPCM_ICP7            0x00000007
#define SPC_STARTDELAY              295110l
#define SPC_BASISTTLTRIG            295120l
#define SPC_TIMEOUT                 295130l
#define SPC_SWL_INFO                295140l
#define SPC_SWD_INFO                295141l
#define SPC_SWD_DOWN                295142l
#define SPC_SWL_EXTRAINFO           295143l
#define SPC_SPECIALCLOCK_ADJUST0    295150l
#define SPC_SPECIALCLOCK_ADJUST1    295151l
#define SPC_SPECIALCLOCK_ADJUST2    295152l
#define SPC_SPECIALCLOCK_ADJUST3    295153l
#define    SPCM_SPECIALCLOCK_ADJUST_SHIFT 1000000
#define SPC_REGACC_CONTMEM          299000l
#define SPC_REGACC_MEMORYUSAGE      299001l
#define SPC_REINITLOGSETTINGS       299998l
#define SPC_LOGDLLCALLS             299999l






// ----- PCK400 -----
#define SPC_FREQUENCE               300000l
#define SPC_DELTAFREQUENCE          300010l
#define SPC_PINHIGH                 300100l
#define SPC_PINLOW                  300110l
#define SPC_PINDELTA                300120l
#define SPC_STOPLEVEL               300200l
#define SPC_PINRELAIS               300210l
#define SPC_EXTERNLEVEL             300300l



// ----- PADCO -----
#define SPC_COUNTER0                310000l
#define SPC_COUNTER1                310001l
#define SPC_COUNTER2                310002l
#define SPC_COUNTER3                310003l
#define SPC_COUNTER4                310004l
#define SPC_COUNTER5                310005l
#define SPC_MODE0                   310100l
#define SPC_MODE1                   310101l
#define SPC_MODE2                   310102l
#define SPC_MODE3                   310103l
#define SPC_MODE4                   310104l
#define SPC_MODE5                   310105l
#define     CM_SINGLE                   1l
#define     CM_MULTI                    2l
#define     CM_POSEDGE                  4l
#define     CM_NEGEDGE                  8l
#define     CM_HIGHPULSE                16l
#define     CM_LOWPULSE                 32l



// ----- PAD1616 -----
#define SPC_SEQUENCERESET           320000l
#define SPC_SEQUENCEADD             320010l
#define     SEQ_IR_10000MV              0l
#define     SEQ_IR_5000MV               1l
#define     SEQ_IR_2000MV               2l
#define     SEQ_IR_1000MV               3l
#define     SEQ_IR_500MV                4l
#define     SEQ_CH0                     0l
#define     SEQ_CH1                     8l
#define     SEQ_CH2                     16l
#define     SEQ_CH3                     24l
#define     SEQ_CH4                     32l
#define     SEQ_CH5                     40l
#define     SEQ_CH6                     48l
#define     SEQ_CH7                     56l
#define     SEQ_CH8                     64l
#define     SEQ_CH9                     72l
#define     SEQ_CH10                    80l
#define     SEQ_CH11                    88l
#define     SEQ_CH12                    96l
#define     SEQ_CH13                    104l
#define     SEQ_CH14                    112l
#define     SEQ_CH15                    120l
#define     SEQ_TRIGGER                 128l
#define     SEQ_START                   256l



// ----- Option CA -----
#define SPC_CA_MODE                 330000l
#define     CAMODE_OFF                  0l
#define     CAMODE_CDM                  1l
#define     CAMODE_KW                   2l
#define     CAMODE_OT                   3l
#define     CAMODE_CDMMUL               4l
#define SPC_CA_TRIGDELAY            330010l
#define SPC_CA_CKDIV                330020l
#define SPC_CA_PULS                 330030l
#define SPC_CA_CKMUL                330040l
#define SPC_CA_DREHZAHLFORMAT       330050l
#define     CADREH_4X4                  0l
#define     CADREH_1X16                 1l
#define SPC_CA_KWINVERT             330060l
#define SPC_CA_OUTA                 330100l
#define SPC_CA_OUTB                 330110l
#define     CAOUT_TRISTATE              0l
#define     CAOUT_LOW                   1l
#define     CAOUT_HIGH                  2l
#define     CAOUT_CDM                   3l
#define     CAOUT_OT                    4l
#define     CAOUT_KW                    5l
#define     CAOUT_TRIG                  6l
#define     CAOUT_CLK                   7l
#define     CAOUT_KW60                  8l
#define     CAOUT_KWGAP                 9l
#define     CAOUT_TRDLY                 10l
#define     CAOUT_INVERT                16l


// ----- Option Sequence Mode (output cards) -----
#define SPC_SEQMODE_STEPMEM0        340000l
// ... 
#define SPC_SEQMODE_STEPMEM8191     348191l

// low part of 64 bit entry
#define     SPCSEQ_SEGMENTMASK      0x0000FFFF
#define     SPCSEQ_NEXTSTEPMASK     0xFFFF0000

// high part of 64 bit entry
#define     SPCSEQ_LOOPMASK         0x000FFFFF
#define     SPCSEQ_ENDLOOPALWAYS    0x00000000
#define     SPCSEQ_ENDLOOPONTRIG    0x40000000
#define     SPCSEQ_END              0x80000000

#define SPC_SEQMODE_AVAILMAXSEGMENT 349900l
#define SPC_SEQMODE_AVAILMAXSTEPS   349901l
#define SPC_SEQMODE_AVAILMAXLOOP    349902l
#define SPC_SEQMODE_AVAILFEATURES   349903l

#define SPC_SEQMODE_MAXSEGMENTS     349910l
#define SPC_SEQMODE_WRITESEGMENT    349920l
#define SPC_SEQMODE_STARTSTEP       349930l
#define SPC_SEQMODE_SEGMENTSIZE     349940l

#define SPC_SEQMODE_STATUS          349950l
#define     SEQSTAT_STEPCHANGE          0x80000000l


// ----- netbox registers -----
#define SPC_NETBOX_TYPE             400000l
#define     NETBOX_SERIES_MASK      0xFF000000
#define     NETBOX_FAMILY_MASK      0x00FF0000
#define     NETBOX_SPEED_MASK       0x0000FF00
#define     NETBOX_CHANNEL_MASK     0x000000FF

#define     NETBOX_SERIES_DN2       0x02000000
#define     NETBOX_SERIES_DN6       0x06000000

#define     NETBOX_FAMILY_20        0x00200000
#define     NETBOX_FAMILY_22        0x00220000
#define     NETBOX_FAMILY_44        0x00440000
#define     NETBOX_FAMILY_46        0x00460000
#define     NETBOX_FAMILY_47        0x00470000
#define     NETBOX_FAMILY_48        0x00480000
#define     NETBOX_FAMILY_49        0x00490000
#define     NETBOX_FAMILY_59        0x00590000
#define     NETBOX_FAMILY_60        0x00600000
#define     NETBOX_FAMILY_65        0x00650000
#define     NETBOX_FAMILY_66        0x00660000
#define     NETBOX_FAMILY_8X        0x00800000
#define     NETBOX_FAMILY_80        0x00800000
#define     NETBOX_FAMILY_81        0x00810000
#define     NETBOX_FAMILY_82        0x00820000
#define     NETBOX_FAMILY_83        0x00830000

#define     NETBOX_SPEED_1          0x00000100
#define     NETBOX_SPEED_2          0x00000200
#define     NETBOX_SPEED_3          0x00000300
#define     NETBOX_SPEED_4          0x00000400
#define     NETBOX_SPEED_5          0x00000500
#define     NETBOX_SPEED_6          0x00000600
#define     NETBOX_SPEED_7          0x00000700
#define     NETBOX_SPEED_8          0x00000800

#define     NETBOX_CHANNELS_2       0x00000002
#define     NETBOX_CHANNELS_4       0x00000004
#define     NETBOX_CHANNELS_6       0x00000006
#define     NETBOX_CHANNELS_8       0x00000008
#define     NETBOX_CHANNELS_10      0x0000000A
#define     NETBOX_CHANNELS_12      0x0000000C
#define     NETBOX_CHANNELS_16      0x00000010
#define     NETBOX_CHANNELS_20      0x00000014
#define     NETBOX_CHANNELS_24      0x00000018
#define     NETBOX_CHANNELS_32      0x00000020
#define     NETBOX_CHANNELS_40      0x00000028
#define     NETBOX_CHANNELS_48      0x00000030

#define SPC_NETBOX_SERIALNO         400001l
#define SPC_NETBOX_PRODUCTIONDATE   400002l
#define SPC_NETBOX_HWVERSION        400003l
#define SPC_NETBOX_SWVERSION        400004l

#define SPC_NETBOX_FEATURES         400005l
#define     NETBOX_FEAT_DCPOWER         0x1
#define     NETBOX_FEAT_BOOTATPOWERON   0x2
#define     NETBOX_FEAT_EMBEDDEDSERVER  0x4

#define SPC_NETBOX_CUSTOM           400006l

#define SPC_NETBOX_WAKEONLAN        400007l
#define SPC_NETBOX_MACADDRESS       400008l
#define SPC_NETBOX_LANIDFLASH       400009l
#define SPC_NETBOX_TEMPERATURE      400010l
#define SPC_NETBOX_SHUTDOWN         400011l
#define SPC_NETBOX_RESTART          400012l
#define SPC_NETBOX_FANSPEED0        400013l
#define SPC_NETBOX_FANSPEED1        400014l
#define SPC_NETBOX_TEMPERATURE_K    400010l // same SPC_NETBOX_TEMPERATURE
#define SPC_NETBOX_TEMPERATURE_C    400015l
#define SPC_NETBOX_TEMPERATURE_F    400016l
#define SPC_NETBOX_TEMPERATURE1_K   400017l
#define SPC_NETBOX_TEMPERATURE1_C   400018l
#define SPC_NETBOX_TEMPERATURE1_F   400019l
#define SPC_NETBOX_TEMPERATURE2_K   400020l
#define SPC_NETBOX_TEMPERATURE2_C   400021l
#define SPC_NETBOX_TEMPERATURE2_F   400022l

// ----- hardware monitor registers -----
#define SPC_MON_V_PCIE_BUS          500000l
#define SPC_MON_V_CONNECTOR         500001l
#define SPC_MON_CARD_PWRSOURCE      500002l
#define     CARD_PWRSOURCE_BUS          0l
#define     CARD_PWRSOURCE_CONNECTOR    1l
#define SPC_MON_V_CARD_IN           500003l
#define SPC_MON_I_CARD_IN           500004l
#define SPC_MON_P_CARD_IN           500005l
#define SPC_MON_V_3V3               500006l
#define SPC_MON_V_2V5               500007l
#define SPC_MON_V_CORE              500008l
#define SPC_MON_V_AVTT              500009l
#define SPC_MON_V_AVCC              500010l
#define SPC_MON_V_MEMVCC            500011l
#define SPC_MON_V_MEMVTT            500012l
#define SPC_MON_V_CP_POS            500013l
#define SPC_MON_V_CP_NEG            500014l

#define SPC_MON_V_5VA               500015l
#define SPC_MON_V_ADCA              500016l
#define SPC_MON_V_ADCD              500017l
#define SPC_MON_V_OP_POS            500018l
#define SPC_MON_V_OP_NEG            500019l
#define SPC_MON_V_COMP_NEG          500020l
#define SPC_MON_V_COMP_POS          500021l

// legacy temperature registers (Kelvin)
#define SPC_MON_T_BASE_CTRL         500022l
#define SPC_MON_T_MODULE_0          500023l
#define SPC_MON_T_MODULE_1          500024l

// new temperature registers for Kelvin (TK), Celsius (TC) or Fahrenheit (TF)
#define SPC_MON_TK_BASE_CTRL         500022l
#define SPC_MON_TK_MODULE_0          500023l
#define SPC_MON_TK_MODULE_1          500024l

#define SPC_MON_TC_BASE_CTRL         500025l
#define SPC_MON_TC_MODULE_0          500026l
#define SPC_MON_TC_MODULE_1          500027l

#define SPC_MON_TF_BASE_CTRL         500028l
#define SPC_MON_TF_MODULE_0          500029l
#define SPC_MON_TF_MODULE_1          500030l

// some more voltages (used on M2p)
#define SPC_MON_V_1V8_BASE           500031l
#define SPC_MON_V_1V8_MOD            500032l
#define SPC_MON_V_MODA_0             500033l
#define SPC_MON_V_MODA_1             500034l
#define SPC_MON_V_MODB_0             500035l
#define SPC_MON_V_MODB_1             500037l

// some more voltages and temperatures (used on M2p.65xx-hv)
#define SPC_MON_TK_MODA_0           500023l // same as SPC_MON_TK_MODULE_0
#define SPC_MON_TK_MODA_1           500038l
#define SPC_MON_TK_MODA_2           500039l
#define SPC_MON_TK_MODA_3           500040l
#define SPC_MON_TK_MODA_4           500041l
#define SPC_MON_TK_MODB_0           500024l // same as SPC_MON_TK_MODULE_1
#define SPC_MON_TK_MODB_1           500042l
#define SPC_MON_TK_MODB_2           500043l
#define SPC_MON_TK_MODB_3           500044l
#define SPC_MON_TK_MODB_4           500045l

#define SPC_MON_TC_MODA_0           500026l // same as SPC_MON_TC_MODULE_0
#define SPC_MON_TC_MODA_1           500046l
#define SPC_MON_TC_MODA_2           500047l
#define SPC_MON_TC_MODA_3           500048l
#define SPC_MON_TC_MODA_4           500049l
#define SPC_MON_TC_MODB_0           500027l // same as SPC_MON_TC_MODULE_1
#define SPC_MON_TC_MODB_1           500050l
#define SPC_MON_TC_MODB_2           500051l
#define SPC_MON_TC_MODB_3           500052l
#define SPC_MON_TC_MODB_4           500053l

#define SPC_MON_TF_MODA_0           500029l // same as SPC_MON_TF_MODULE_0
#define SPC_MON_TF_MODA_1           500054l
#define SPC_MON_TF_MODA_2           500055l
#define SPC_MON_TF_MODA_3           500056l
#define SPC_MON_TF_MODA_4           500057l
#define SPC_MON_TF_MODB_0           500030l // same as SPC_MON_TF_MODULE_1
#define SPC_MON_TF_MODB_1           500058l
#define SPC_MON_TF_MODB_2           500059l
#define SPC_MON_TF_MODB_3           500060l
#define SPC_MON_TF_MODB_4           500061l

#define SPC_MON_I_MODA_0            500062l
#define SPC_MON_I_MODA_1            500063l
#define SPC_MON_I_MODA_2            500064l
#define SPC_MON_I_MODA_3            500065l
#define SPC_MON_I_MODB_0            500066l
#define SPC_MON_I_MODB_1            500067l
#define SPC_MON_I_MODB_2            500068l
#define SPC_MON_I_MODB_3            500069l

#define SPC_MON_MOD_FAULT           500070l
#define SPC_CLR_MOD_FAULT           500071l

// power section temperature registers for Kelvin (TK), Celsius (TC) or Fahrenheit (TF)
#define SPC_MON_TK_MODA_5           500072l
#define SPC_MON_TK_MODB_5           500073l

#define SPC_MON_TC_MODA_5           500074l
#define SPC_MON_TC_MODB_5           500075l

#define SPC_MON_TF_MODA_5           500076l
#define SPC_MON_TF_MODB_5           500077l

// mask with available monitor registers
#define SPC_AVAILMONITORS            510000l
#define     SPCM_MON_T_BASE_CTRL        0x0000000000000001ULL
#define     SPCM_MON_T_MODULE_0         0x0000000000000002ULL
#define     SPCM_MON_T_MODULE_1         0x0000000000000004ULL

#define     SPCM_MON_V_PCIE_BUS         0x0000000000000010ULL
#define     SPCM_MON_V_CONNECTOR        0x0000000000000020ULL
#define     SPCM_MON_CARD_PWRSOURCE     0x0000000000000040ULL
#define     SPCM_MON_V_CARD_IN          0x0000000000000080ULL
#define     SPCM_MON_I_CARD_IN          0x0000000000000100ULL
#define     SPCM_MON_P_CARD_IN          0x0000000000000200ULL
#define     SPCM_MON_V_3V3              0x0000000000000400ULL
#define     SPCM_MON_V_2V5              0x0000000000000800ULL
#define     SPCM_MON_V_CORE             0x0000000000001000ULL
#define     SPCM_MON_V_AVTT             0x0000000000002000ULL
#define     SPCM_MON_V_AVCC             0x0000000000004000ULL
#define     SPCM_MON_V_MEMVCC           0x0000000000008000ULL
#define     SPCM_MON_V_MEMVTT           0x0000000000010000ULL
#define     SPCM_MON_V_CP_POS           0x0000000000020000ULL
#define     SPCM_MON_V_CP_NEG           0x0000000000040000ULL
#define     SPCM_MON_V_5VA              0x0000000000080000ULL
#define     SPCM_MON_V_ADCA             0x0000000000100000ULL
#define     SPCM_MON_V_ADCD             0x0000000000200000ULL
#define     SPCM_MON_V_OP_POS           0x0000000000400000ULL
#define     SPCM_MON_V_OP_NEG           0x0000000000800000ULL
#define     SPCM_MON_V_COMP_NEG         0x0000000001000000ULL
#define     SPCM_MON_V_COMP_POS         0x0000000002000000ULL
#define     SPCM_MON_V_1V8_BASE         0x0000000004000000ULL
#define     SPCM_MON_V_1V8_MOD          0x0000000008000000ULL

#define     SPCM_MON_V_MODA_0           0x0000000010000000ULL
#define     SPCM_MON_V_MODA_1           0x0000000020000000ULL
#define     SPCM_MON_V_MODB_0           0x0000000040000000ULL
#define     SPCM_MON_V_MODB_1           0x0000000080000000ULL

#define     SPCM_MON_T_MODA_0           0x0000000000000002ULL // same as SPCM_MON_T_MODULE_0
#define     SPCM_MON_T_MODA_1           0x0000000100000000ULL
#define     SPCM_MON_T_MODA_2           0x0000000200000000ULL
#define     SPCM_MON_T_MODA_3           0x0000000400000000ULL
#define     SPCM_MON_T_MODA_4           0x0000000800000000ULL
#define     SPCM_MON_T_MODB_0           0x0000000000000004ULL // same as SPCM_MON_T_MODULE_1
#define     SPCM_MON_T_MODB_1           0x0000001000000000ULL
#define     SPCM_MON_T_MODB_2           0x0000002000000000ULL
#define     SPCM_MON_T_MODB_3           0x0000004000000000ULL
#define     SPCM_MON_T_MODB_4           0x0000008000000000ULL

#define     SPCM_MON_I_MODA_0           0x0000010000000000ULL
#define     SPCM_MON_I_MODA_1           0x0000020000000000ULL
#define     SPCM_MON_I_MODA_2           0x0000040000000000ULL
#define     SPCM_MON_I_MODA_3           0x0000080000000000ULL
#define     SPCM_MON_I_MODB_0           0x0000100000000000ULL
#define     SPCM_MON_I_MODB_1           0x0000200000000000ULL
#define     SPCM_MON_I_MODB_2           0x0000300000000000ULL
#define     SPCM_MON_I_MODB_3           0x0000400000000000ULL

#define     SPCM_MON_T_MODA_5           0x0000800000000000ULL
#define     SPCM_MON_T_MODB_5           0x0001000000000000ULL


// ----- re-located multi-purpose i/o related registers -----
#define SPC_X0_READFEATURES         600000l
#define SPC_X1_READFEATURES         600001l
#define SPC_X2_READFEATURES         600002l
#define SPC_X3_READFEATURES         600003l
#define SPC_X4_READFEATURES         600004l
#define SPC_X5_READFEATURES         600005l
#define SPC_X6_READFEATURES         600006l
#define SPC_X7_READFEATURES         600007l
#define SPC_X8_READFEATURES         600008l
#define SPC_X9_READFEATURES         600009l
#define SPC_X10_READFEATURES        600010l
#define SPC_X11_READFEATURES        600011l
#define SPC_X12_READFEATURES        600012l
#define SPC_X13_READFEATURES        600013l
#define SPC_X14_READFEATURES        600014l
#define SPC_X15_READFEATURES        600015l
#define SPC_X16_READFEATURES        600016l
#define SPC_X17_READFEATURES        600017l
#define SPC_X18_READFEATURES        600018l
#define SPC_X19_READFEATURES        600019l
#define     SPCM_XFEAT_TERM             0x00000001
#define     SPCM_XFEAT_HIGHIMP          0x00000002
#define     SPCM_XFEAT_DCCOUPLING       0x00000004
#define     SPCM_XFEAT_ACCOUPLING       0x00000008
#define     SPCM_XFEAT_SE               0x00000010
#define     SPCM_XFEAT_DIFF             0x00000020
#define     SPCM_XFEAT_PROGTHRESHOLD    0x00000040

#define SPC_X0_TERM                600100l
#define SPC_X1_TERM                600101l
#define SPC_X2_TERM                600102l
#define SPC_X3_TERM                600103l
#define SPC_X4_TERM                600104l
#define SPC_X5_TERM                600105l
#define SPC_X6_TERM                600106l
#define SPC_X7_TERM                600107l
#define SPC_X8_TERM                600108l
#define SPC_X9_TERM                600109l
#define SPC_X10_TERM               600110l
#define SPC_X11_TERM               600111l
#define SPC_X12_TERM               600112l
#define SPC_X13_TERM               600113l
#define SPC_X14_TERM               600114l
#define SPC_X15_TERM               600115l
#define SPC_X16_TERM               600116l
#define SPC_X17_TERM               600117l
#define SPC_X18_TERM               600118l
#define SPC_X19_TERM               600119l

#define SPCM_X0_MODE                600200l
#define SPCM_X1_MODE                600201l
#define SPCM_X2_MODE                600202l
#define SPCM_X3_MODE                600203l
#define SPCM_X4_MODE                600204l
#define SPCM_X5_MODE                600205l
#define SPCM_X6_MODE                600206l
#define SPCM_X7_MODE                600207l
#define SPCM_X8_MODE                600208l
#define SPCM_X9_MODE                600209l
#define SPCM_X10_MODE               600210l
#define SPCM_X11_MODE               600211l
#define SPCM_X12_MODE               600212l
#define SPCM_X13_MODE               600213l
#define SPCM_X14_MODE               600214l
#define SPCM_X15_MODE               600215l
#define SPCM_X16_MODE               600216l
#define SPCM_X17_MODE               600217l
#define SPCM_X18_MODE               600218l
#define SPCM_X19_MODE               600219l

#define SPCM_X0_AVAILMODES          600300l
#define SPCM_X1_AVAILMODES          600301l
#define SPCM_X2_AVAILMODES          600302l
#define SPCM_X3_AVAILMODES          600303l
#define SPCM_X4_AVAILMODES          600304l
#define SPCM_X5_AVAILMODES          600305l
#define SPCM_X6_AVAILMODES          600306l
#define SPCM_X7_AVAILMODES          600307l
#define SPCM_X8_AVAILMODES          600308l
#define SPCM_X9_AVAILMODES          600309l
#define SPCM_X10_AVAILMODES         600310l
#define SPCM_X11_AVAILMODES         600311l
#define SPCM_X12_AVAILMODES         600312l
#define SPCM_X13_AVAILMODES         600313l
#define SPCM_X14_AVAILMODES         600314l
#define SPCM_X15_AVAILMODES         600315l
#define SPCM_X16_AVAILMODES         600316l
#define SPCM_X17_AVAILMODES         600317l
#define SPCM_X18_AVAILMODES         600318l
#define SPCM_X19_AVAILMODES         600319l
// for definitions of the available modes see section at SPCM_LEGACY_X0_MODE above


// ----- Hardware registers (debug use only) -----
#define SPC_REG0x00                 900000l
#define SPC_REG0x02                 900010l
#define SPC_REG0x04                 900020l
#define SPC_REG0x06                 900030l
#define SPC_REG0x08                 900040l
#define SPC_REG0x0A                 900050l
#define SPC_REG0x0C                 900060l
#define SPC_REG0x0E                 900070l

#define SPC_DEBUGREG0               900100l
#define SPC_DEBUGREG15              900115l
#define SPC_DEBUGVALUE0             900200l
#define SPC_DEBUGVALUE15            900215l

#define SPC_MI_ISP                  901000l
#define     ISP_TMS_0                   0l
#define     ISP_TMS_1                   1l
#define     ISP_TDO_0                   0l
#define     ISP_TDO_1                   2l


#define SPC_EE_RWAUTH               901100l
#define SPC_EE_REG                  901110l
#define SPC_EE_RESETCOUNTER         901120l

// ----- Test Registers -----
#define SPC_TEST_BASE               902000l
#define SPC_TEST_LOCAL_START        902100l
#define SPC_TEST_LOCAL_END          902356l
#define SPC_TEST_PLX_START          902400l
#define SPC_TEST_PLX_END            902656l

// 9012xx not usable
// 901900 not usable
// 903000 not usable
// 91xxxx not usable

// ----- used by GetErrorInfo to mark errors in other functions than SetParam/GetParam -----
#define SPC_FUNCTION_DEFTRANSFER 100000000l
