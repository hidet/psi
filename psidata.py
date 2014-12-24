

pulse_dir=["20141026_A",   "20141026_B",   "20141026_C",   "20141027_A",   "20141027_B",
           "20141027_C",   "20141027_D",   "20141028_A",   "20141028_B",   "20141028_C",
           "20141028_D",   "20141028_E",   "20141028_F",   "20141028_G",   "20141028_H",
           "20141028_I",   "20141028_J",   "20141029_A",   "20141029_B",   "20141029_C",
           "20141029_D",   "20141029_E",   "20141029_F",   "20141029_G",   "20141029_H",
           "20141029_I",   "20141029_J",   "20141029_K",   "20141029_L",   "20141029_M",
           "20141029_N",   "20141029_O",   "20141029_P",   "20141029_Q",   "20141029_R",
           "20141029_S",   "20141029_T",   "20141029_U",   "20141029_V",   "20141029_W",
           "20141029_X",   "20141029_Y",   "20141029_Z",   "20141029_AA",  "20141029_AB",
           "20141029_AC",  "20141029_AD",  "20141029_AE",  "20141029_AF",  "20141029_AG",
           "20141029_AH",  "20141030_A",   "20141030_B",   "20141030_C",   "20141030_D",
           "20141030_E",   "20141030_F",   "20141030_G",   "20141030_H",   "20141030_I",
           "20141030_J",   "20141030_K",   "20141030_L",   "20141030_M",   "20141031_A",
           "20141031_B",   "20141031_C",   "20141031_D",   "20141031_E",   "20141031_F",
           "20141031_G",   "20141031_H",   "20141031_I",   "20141031_J",   "20141031_K",
           "20141031_L",   "20141031_M",   "20141031_N",   "20141031_O",   "20141031_P",
           "20141031_Q",   "20141031_R",   "20141031_S",   "20141101_A",   "20141101_B",
           "20141101_C",   "20141101_D",   "20141101_E",   "20141101_F",   "20141101_G",
           "20141101_H",   "20141101_I",   "20141101_J",   "20141101_K",   "20141101_L",
           "20141101_M",   "20141101_N",   "20141101_O",   "20141101_P",   "20141101_Q",
           "20141101_R",   "20141101_S",   "20141101_T",   "20141101_U",   "20141101_V",
           "20141102_A",   "20141102_B",   "20141102_C",   "20141102_D",   "20141102_E",
           "20141102_F",   "20141102_G",   "20141102_H",   "20141102_I",   "20141102_J",
           "20141102_K",   "20141102_L",   "20141102_M",   "20141102_N",   "20141102_O",
           "20141102_P",   "20141102_Q",   "20141102_R",   "20141102_S",   "20141102_T",
           "20141102_U",   "20141102_V",   "20141102_W",   "20141102_X",   "20141103_A",
           "20141103_B",   "20141103_C",   "20141103_D",   "20141103_E",   "20141103_F",
           "20141103_G",   "20141103_H",   "20141103_I",   "20141103_J",   "20141103_K",
           "20141103_L",   "20141103_M",   "20141103_N",   "20141104_A",   "20141104_B",
           "20141104_C",   "20141104_D",   "20141104_E",   "20141104_F",   "20141104_G",
           "20141104_H",   "20141104_I",   "20141104_J",   "20141104_K",   "20141104_L",
           "20141104_M",   "20141104_N",   "20141105_A",   "20141105_B",   "20141105_C",
           "20141105_D",   "20141105_E",   "20141105_F",   "20141105_G",   "20141105_H",
           "20141105_I",   "20141105_J",   "20141105_K",   "20141105_L",   "20141105_M",
           "20141105_O",   "20141105_P",   "20141105_Q",   "20141105_R",   "20141105_S"]



noise_name = ["A",   "A",   "A",   "A",   "A",   "A",   "A",   "A",  "A",   "A",
              "A",   "A",   "A",   "G",   "G",   "G",   "G",   "A",   "A",   "A",
              "A",   "A",   "F",   "F",   "F",   "F",   "F",   "F",   "F",   "F",
              "F",   "F",   "P",   "P",   "P",   "P",   "P",   "P",   "P",   "P",
              "P",   "P",   "P",   "P",   "P",   "P",   "P",   "AE",  "AE",  "AE",
              "AE",  "A",   "A",   "A",   "A",   "A",   "F",   "F",   "F",   "F",
              "F",   "F",   "F",   "F",   "F",   "F",   "F",   "F",   "F",   "F",
              "F",   "F",   "F",   "F",   "F",   "F",   "F",   "F",   "F",   "F",
              "F",   "F",   "F",   "A",   "A",   "A",   "A",   "A",   "A",   "A",
              "H",   "I",   "I",   "I",   "L",   "L",   "L",   "L",   "L",   "L",
              "L",   "L",   "L",   "L",   "L",   "L",   "L",   "L",   "L",   "L",
              "L",   "L",   "L",   "L",   "L",   "L",   "L",   "L",   "L",   "L",
              "L",   "L",   "R",   "R",   "R",   "R",   "R",   "R",   "R",   "R",
              "R",   "R",   "R",   "R",   "R",   "R",   "R",   "R",   "R",   "R",
              "R",   "R",   "R",   "R",   "R",   "C",   "C",   "C",   "C",   "C",
              "C",   "C",   "C",   "C",   "C",   "C",   "C",   "C",   "C",   "C",
              "C",   "C",   "C",   "C",   "C",   "C",   "C",   "C",   "C",   "C",
              "S",   "S",   "S",   "S",   "S"]



noise_date = ["20141026",   "20141026",   "20141026",   "20141027",   "20141027",
              "20141027",   "20141027",   "20141028",   "20141028",   "20141028",
              "20141028",   "20141028",   "20141028",   "20141028",   "20141028",
              "20141028",   "20141028",   "20141029",   "20141029",   "20141029",
              "20141029",   "20141029",   "20141029",   "20141029",   "20141029",
              "20141029",   "20141029",   "20141029",   "20141029",   "20141029",
              "20141029",   "20141029",   "20141029",   "20141029",   "20141029",
              "20141029",   "20141029",   "20141029",   "20141029",   "20141029",
              "20141029",   "20141029",   "20141029",   "20141029",   "20141029",
              "20141029",   "20141029",   "20141029",   "20141029",   "20141029",
              "20141029",   "20141030",   "20141030",   "20141030",   "20141030",
              "20141030",   "20141030",   "20141030",   "20141030",   "20141030",
              "20141030",   "20141030",   "20141030",   "20141030",   "20141030",
              "20141030",   "20141030",   "20141030",   "20141030",   "20141030",
              "20141030",   "20141030",   "20141030",   "20141030",   "20141030",
              "20141030",   "20141030",   "20141030",   "20141030",   "20141030",
              "20141030",   "20141030",   "20141030",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141101",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141101",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141101",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141101",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141101",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141101",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141101",   "20141101",   "20141101",
              "20141101",   "20141101",   "20141102",   "20141102",   "20141102",
              "20141102",   "20141102",   "20141102",   "20141102",   "20141103",
              "20141103",   "20141103",   "20141103",   "20141103",   "20141103",
              "20141103",   "20141103",   "20141103",   "20141103",   "20141103",
              "20141103",   "20141103",   "20141103",   "20141103",   "20141103",
              "20141104",   "20141104",   "20141104",   "20141104",   "20141104",
              "20141104",   "20141104",   "20141104",   "20141104",   "20141104",
              "20141104",   "20141104",   "20141104",   "20141104",   "20141104",
              "20141104",   "20141104",   "20141104",   "20141104",   "20141104",
              "20141104",   "20141104",   "20141104",   "20141104",   "20141104",
              "20141105",   "20141105",   "20141105",   "20141105",   "20141105"]   


ana_list = ["",       "Mn",     "RbBr",   "",       "lowE",   "",       "lowE",   "",       "lowE",   "Se",
           "Se",     "RbBr",   "Mn",     "",       "Mn",     "GaAs",   "highE",  "",       "lowE",   "lowE",
           "Off",       "lowE",   "",       "lowE",   "lowE",   "lowE",   "lowE",   "",       "lowE",   "Mn",
           "Mn",     "Mn",     "",       "Mn",     "Mn",     "Mn",     "Mn",     "Off",    "Mn",     "Mn",
           "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "",       "",       "Mn",     "Mn",
           "Mn",     "",       "lowE",   "lowE",   "lowE",   "lowE",   "",       "Mn",     "lowE",   "lowE",
           "lowE",   "lowE",   "lowE",   "lowE",   "lowE",   "Off",   "lowE",   "Off",   "lowE",   "Off",
           "lowE",   "Off",   "CrCo",   "CrCo",   "CrCo",   "CrCo",   "CrCo",   "CrCo",   "CrCo",   "CrCo",
           "CrCo",   "CrCo",   "CrCo",   "",       "CrCo",   "CrCo",   "CrCo",   "CrCo",   "CrCo",   "CrCo",
           "",       "",       "CrCo",   "CrCo",   "",       "CrCo",   "CrCo",   "CrCo",   "CrCo",   "CrCo",
           "lowE",   "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",
           "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "Mn",     "CrCo",   "CrCo",
           "CrCo",   "CrCo",   "",       "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",
           "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",
           "highE",  "highE",  "highE",  "highE",  "highE",  "",       "highE",  "highE",  "highE",  "highE",
           "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",
           "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",  "highE",
           "Mn",     "Mn",     "Mn",     "Mn",     ""]



def getnoiseana(dir_pl):
    ind = pulse_dir.index(dir_pl)
    ndate = noise_date[ind]
    nname = noise_name[ind]
    ana = ana_list[ind]
    return ndate, nname, ana



'''
"20141026_A"  "20141026"     "A"      ""       
"20141026_B"  "20141026"     "A"      "Mn"     
"20141026_C"  "20141026"     "A"      "RbBr"   
"20141027_A"  "20141027"     "A"      ""       
"20141027_B"  "20141027"     "A"      "lowE"   
"20141027_C"  "20141027"     "A"      ""       
"20141027_D"  "20141027"     "A"      "lowE"   
"20141028_A"  "20141028"     "A"      ""       
"20141028_B"  "20141028"     "A"      "lowE"   
"20141028_C"  "20141028"     "A"      "Se"     
"20141028_D"  "20141028"     "A"      "Se"     
"20141028_E"  "20141028"     "A"      "RbBr"   
"20141028_F"  "20141028"     "A"      "Mn"     
"20141028_G"  "20141028"     "G"      ""       
"20141028_H"  "20141028"     "G"      "Mn"     
"20141028_I"  "20141028"     "G"      "GaAs"   
"20141028_J"  "20141028"     "G"      "highE"  
"20141029_A"  "20141029"     "A"      ""       
"20141029_B"  "20141029"     "A"      "lowE"   
"20141029_C"  "20141029"     "A"      "lowE"   
"20141029_D"  "20141029"     "A"      ""       
"20141029_E"  "20141029"     "A"      "lowE"   
"20141029_F"  "20141029"     "F"      ""       
"20141029_G"  "20141029"     "F"      "lowE"   
"20141029_H"  "20141029"     "F"      "lowE"   
"20141029_I"  "20141029"     "F"      "lowE"   
"20141029_J"  "20141029"     "F"      "lowE"   
"20141029_K"  "20141029"     "F"      ""       
"20141029_L"  "20141029"     "F"      "lowE"   
"20141029_M"  "20141029"     "F"      "Mn"     
"20141029_N"  "20141029"     "F"      "Mn"     
"20141029_O"  "20141029"     "F"      "Mn"     
"20141029_P"  "20141029"     "P"      ""       
"20141029_Q"  "20141029"     "P"      "Mn"     
"20141029_R"  "20141029"     "P"      "Mn"     
"20141029_S"  "20141029"     "P"      "Mn"     
"20141029_T"  "20141029"     "P"      "Mn"     
"20141029_U"  "20141029"     "P"      "Mn"     
"20141029_V"  "20141029"     "P"      "Mn"     
"20141029_W"  "20141029"     "P"      "Mn"     
"20141029_X"  "20141029"     "P"      "Mn"     
"20141029_Y"  "20141029"     "P"      "Mn"     
"20141029_Z"  "20141029"     "P"      "Mn"     
"20141029_AA" "20141029"     "P"      "Mn"     
"20141029_AB" "20141029"     "P"      "Mn"     
"20141029_AC" "20141029"     "P"      "Mn"     
"20141029_AD" "20141029"     "P"      ""       
"20141029_AE" "20141029"     "AE"     ""       
"20141029_AF" "20141029"     "AE"     "Mn"     
"20141029_AG" "20141029"     "AE"     "Mn"     
"20141029_AH" "20141029"     "AE"     "Mn"     
"20141030_A"  "20141030"     "A"      ""       
"20141030_B"  "20141030"     "A"      "lowE"   
"20141030_C"  "20141030"     "A"      "lowE"   
"20141030_D"  "20141030"     "A"      "lowE"   
"20141030_E"  "20141030"     "A"      "lowE"   
"20141030_F"  "20141030"     "F"      ""       
"20141030_G"  "20141030"     "F"      "Mn"     
"20141030_H"  "20141030"     "F"      "lowE"   
"20141030_I"  "20141030"     "F"      "lowE"   
"20141030_J"  "20141030"     "F"      "lowE"   
"20141030_K"  "20141030"     "F"      "lowE"   
"20141030_L"  "20141030"     "F"      "lowE"   
"20141030_M"  "20141030"     "F"      "lowE"   
"20141031_A"  "20141030"     "F"      "lowE"   
"20141031_B"  "20141030"     "F"      "lowE"   
"20141031_C"  "20141030"     "F"      "lowE"   
"20141031_D"  "20141030"     "F"      "lowE"   
"20141031_E"  "20141030"     "F"      "lowE"   
"20141031_F"  "20141030"     "F"      "lowE"   
"20141031_G"  "20141030"     "F"      "lowE"   
"20141031_H"  "20141030"     "F"      "lowE"   
"20141031_I"  "20141030"     "F"      "CrCo"   
"20141031_J"  "20141030"     "F"      "CrCo"   
"20141031_K"  "20141030"     "F"      "CrCo"   
"20141031_L"  "20141030"     "F"      "CrCo"   
"20141031_M"  "20141030"     "F"      "CrCo"   
"20141031_N"  "20141030"     "F"      "CrCo"   
"20141031_O"  "20141030"     "F"      "CrCo"   
"20141031_P"  "20141030"     "F"      "CrCo"   
"20141031_Q"  "20141030"     "F"      "CrCo"   
"20141031_R"  "20141030"     "F"      "CrCo"   
"20141031_S"  "20141030"     "F"      "CrCo"   
"20141101_A"  "20141101"     "A"      ""       
"20141101_B"  "20141101"     "A"      "CrCo"   
"20141101_C"  "20141101"     "A"      "CrCo"   
"20141101_D"  "20141101"     "A"      "CrCo"   
"20141101_E"  "20141101"     "A"      "CrCo"   
"20141101_F"  "20141101"     "A"      "CrCo"   
"20141101_G"  "20141101"     "A"      "CrCo"   
"20141101_H"  "20141101"     "H"      ""       
"20141101_I"  "20141101"     "I"      ""       
"20141101_J"  "20141101"     "I"      "CrCo"   
"20141101_K"  "20141101"     "I"      "CrCo"   
"20141101_L"  "20141101"     "L"      ""       
"20141101_M"  "20141101"     "L"      "CrCo"   
"20141101_N"  "20141101"     "L"      "CrCo"   
"20141101_O"  "20141101"     "L"      "CrCo"   
"20141101_P"  "20141101"     "L"      "CrCo"   
"20141101_Q"  "20141101"     "L"      "CrCo"   
"20141101_R"  "20141101"     "L"      "lowE"   
"20141101_S"  "20141101"     "L"      "Mn"     
"20141101_T"  "20141101"     "L"      "Mn"     
"20141101_U"  "20141101"     "L"      "Mn"     
"20141101_V"  "20141101"     "L"      "Mn"     
"20141102_A"  "20141101"     "L"      "Mn"     
"20141102_B"  "20141101"     "L"      "Mn"     
"20141102_C"  "20141101"     "L"      "Mn"     
"20141102_D"  "20141101"     "L"      "Mn"     
"20141102_E"  "20141101"     "L"      "Mn"     
"20141102_F"  "20141101"     "L"      "Mn"     
"20141102_G"  "20141101"     "L"      "Mn"     
"20141102_H"  "20141101"     "L"      "Mn"     
"20141102_I"  "20141101"     "L"      "Mn"     
"20141102_J"  "20141101"     "L"      "Mn"     
"20141102_K"  "20141101"     "L"      "Mn"     
"20141102_L"  "20141101"     "L"      "Mn"     
"20141102_M"  "20141101"     "L"      "Mn"     
"20141102_N"  "20141101"     "L"      "CrCo"   
"20141102_O"  "20141101"     "L"      "CrCo"   
"20141102_P"  "20141101"     "L"      "CrCo"   
"20141102_Q"  "20141101"     "L"      "CrCo"   
"20141102_R"  "20141102"     "R"      ""       
"20141102_S"  "20141102"     "R"      "highE"  
"20141102_T"  "20141102"     "R"      "highE"  
"20141102_U"  "20141102"     "R"      "highE"  
"20141102_V"  "20141102"     "R"      "highE"  
"20141102_W"  "20141102"     "R"      "highE"  
"20141102_X"  "20141102"     "R"      "highE"  
"20141103_A"  "20141103"     "R"      "highE"  
"20141103_B"  "20141103"     "R"      "highE"  
"20141103_C"  "20141103"     "R"      "highE"  
"20141103_D"  "20141103"     "R"      "highE"  
"20141103_E"  "20141103"     "R"      "highE"  
"20141103_F"  "20141103"     "R"      "highE"  
"20141103_G"  "20141103"     "R"      "highE"  
"20141103_H"  "20141103"     "R"      "highE"  
"20141103_I"  "20141103"     "R"      "highE"  
"20141103_J"  "20141103"     "R"      "highE"  
"20141103_K"  "20141103"     "R"      "highE"  
"20141103_L"  "20141103"     "R"      "highE"  
"20141103_M"  "20141103"     "R"      "highE"  
"20141103_N"  "20141103"     "R"      "highE"  
"20141104_A"  "20141103"     "R"      "highE"  
"20141104_B"  "20141103"     "R"      "highE"  
"20141104_C"  "20141104"     "C"      ""       
"20141104_D"  "20141104"     "C"      "highE"  
"20141104_E"  "20141104"     "C"      "highE"  
"20141104_F"  "20141104"     "C"      "highE"  
"20141104_G"  "20141104"     "C"      "highE"  
"20141104_H"  "20141104"     "C"      "highE"  
"20141104_I"  "20141104"     "C"      "highE"  
"20141104_J"  "20141104"     "C"      "highE"  
"20141104_K"  "20141104"     "C"      "highE"  
"20141104_L"  "20141104"     "C"      "highE"  
"20141104_M"  "20141104"     "C"      "highE"  
"20141104_N"  "20141104"     "C"      "highE"  
"20141105_A"  "20141104"     "C"      "highE"  
"20141105_B"  "20141104"     "C"      "highE"  
"20141105_C"  "20141104"     "C"      "highE"  
"20141105_D"  "20141104"     "C"      "highE"  
"20141105_E"  "20141104"     "C"      "highE"  
"20141105_F"  "20141104"     "C"      "highE"  
"20141105_G"  "20141104"     "C"      "highE"  
"20141105_H"  "20141104"     "C"      "highE"  
"20141105_I"  "20141104"     "C"      "highE"  
"20141105_J"  "20141104"     "C"      "highE"  
"20141105_K"  "20141104"     "C"      "highE"  
"20141105_L"  "20141104"     "C"      "highE"  
"20141105_M"  "20141104"     "C"      "highE"  
"20141105_O"  "20141105"     "S"      "Mn"     
"20141105_P"  "20141105"     "S"      "Mn"     
"20141105_Q"  "20141105"     "S"      "Mn"     
"20141105_R"  "20141105"     "S"      "Mn"     
"20141105_S"  "20141105"     "S"      ""       
'''
