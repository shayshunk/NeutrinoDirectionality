// Period Detector Configs
#include <vector>
using std::vector;

vector<int> excludeList1{0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  18,
                         21,  23,  24,  25,  26,  27,  28,  32,  40,  41,  42,  44,  55,  56,  68,  69,
                         70,  73,  79,  83,  84,  97,  98,  102, 107, 111, 112, 122, 125, 126, 127, 130,
                         139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153};

vector<int> excludeList2{0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  21,  23,  24,  25,  26,
                         27,  28,  32,  34,  40,  41,  42,  44,  52,  55,  56,  69,  70,  79,  83,  84,  86,  97,  98,  111,
                         112, 115, 122, 125, 126, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153};

vector<int> excludeList3{0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  18,  21,
                         23,  24,  25,  26,  27,  28,  31,  32,  34,  40,  41,  42,  44,  48,  52,  55,  56,
                         63,  69,  70,  79,  83,  84,  86,  87,  97,  98,  111, 112, 115, 122, 125, 126, 127,
                         139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153};

vector<int> excludeList4{0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  18,  21,  23,  24,
                         25,  26,  27,  28,  31,  32,  34,  40,  41,  42,  43,  44,  46,  47,  48,  50,  52,  55,  56,
                         63,  69,  70,  73,  79,  83,  84,  86,  87,  97,  98,  111, 112, 115, 122, 125, 126, 127, 128,
                         133, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153};

vector<int> excludeList5{0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  18,  21,  23,  24,  25,
                         26,  27,  28,  29,  31,  32,  34,  36,  40,  41,  42,  43,  44,  46,  47,  48,  50,  52,  55,  56,
                         60,  63,  69,  70,  73,  79,  83,  84,  86,  87,  94,  97,  98,  111, 112, 115, 121, 122, 125, 126,
                         127, 128, 133, 136, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153};

vector<vector<int>> excludeList = {excludeList1, excludeList2, excludeList3, excludeList4, excludeList5};
vector<vector<int>> detectorConfig;