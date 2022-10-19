// constants
#define maxf 8193 // max file name length 

// filename
char *inputfile, *outputfile;

// sequence information
int seq_type = 0; // 1 DNA, 2 Protein, 0 AUTO
int seq_num; // the number of sequences
int seq_maxlen; // the max length of sequences
int maxmemory; // max memory for CD-HIT

// arguments on MAFFT
int fftWinSize, fftthreshold, BLOSUM, ppenalty, ppenalty_dist, alignband, threads, nmax_shift, printdebug;

// arguments on CD-HIT
double cdhitsim;

// center sequence file
char centerfile[maxf];

// command stream
char cmdstr[maxf], cmdstr2[maxf], cmdstr3[maxf];

// tmp directory
char tmpdir[maxf] = "./swap/", *readtmpdir;
int tmpinthisdir;

char orderprotein[20] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
int BLOSUM62[20][20] = {
{  483,  -40,  -52,  -74,   58,   19,   13,  114,  -61,  -31,  -46,   26,    7, -119,   18,  207,   94, -150,  -75,   80},
{  -40,  635,   55,  -59, -234,  194,   87, -128,   74, -195, -113,  305,  -36, -175, -108,   23,  -12, -165,  -68, -147},
{  -52,   55,  653,  223, -163,   98,   72,   57,  155, -218, -233,   81, -113, -195,  -98,  157,   94, -264, -106, -184},
{  -74,  -59,  223,  665, -241,   68,  247,  -31,  -12, -208, -256,   30, -202, -244,  -47,   73,   -5, -315, -203, -210},
{   58, -234, -163, -241,  940, -186, -256, -147, -195,  -22,  -27, -199,  -41, -135, -176,   13,   13, -128, -138,   19},
{   19,  194,   98,   68, -186,  617,  280,  -77,  142, -173, -111,  223,   57, -212,  -27,   88,   32,  -93,  -41, -118},
{   13,   87,   72,  247, -256,  280,  579, -109,   86, -215, -181,  174,  -98, -215,  -11,   84,   14, -180, -100, -141},
{  114, -128,   57,  -31, -147,  -77, -109,  644, -102, -267, -257,  -51, -164, -206, -111,   69,  -56, -146, -200, -210},
{  -61,   74,  155,  -12, -195,  142,   86, -102,  835, -219, -175,   28,  -54,  -23, -114,   12,  -67, -131,  264, -208},
{  -31, -195, -218, -208,  -22, -173, -215, -267, -219,  491,  247, -164,  209,   83, -172, -132,   28, -155,  -33,  348},
{  -46, -113, -233, -256,  -27, -111, -181, -257, -175,  247,  476, -142,  294,  139, -182, -141,  -19,  -62,   -6,  175},
{   26,  305,   81,   30, -199,  223,  174,  -51,   28, -164, -142,  540,  -34, -204,   -1,   78,   33, -191,  -80, -123},
{    7,  -36, -113, -202,  -41,   57,  -98, -164,  -54,  209,  294,  -34,  627,  100, -144,  -47,   33,  -42,    1,  166},
{ -119, -175, -195, -244, -135, -212, -215, -206,  -23,   83,  139, -204,  100,  691, -255, -134, -108,  188,  387,   15},
{   18, -108,  -98,  -47, -176,  -27,  -11, -111, -114, -172, -182,   -1, -144, -255,  821,   19,   -7, -260, -188, -132},
{  207,   23,  157,   73,   13,   88,   84,   69,   12, -132, -141,   78,  -47, -134,   19,  479,  234, -172,  -67,  -63},
{   94,  -12,   94,   -5,   13,   32,   14,  -56,  -67,   28,  -19,   33,   33, -108,   -7,  234,  544, -140,  -59,   93},
{ -150, -165, -264, -315, -128,  -93, -180, -146, -131, -155,  -62, -191,  -42,  188, -260, -172, -140, 1129,  309, -180},
{  -75,  -68, -106, -203, -138,  -41, -100, -200,  264,  -33,   -6,  -80,    1,  387, -188,  -67,  -59,  309,  745,  -20},
{   80, -147, -184, -210,   19, -118, -141, -210, -208,  348,  175, -123,  166,   15, -132,  -63,   93, -180,  -20,  468}
};

char orderDNA[4] = {'a', 't', 'c', 'g'};
int trans[4][4] = {{600, 128, -364, -364}, {128, 600, -364, -364}, {-364, -364, 600, 128}, {-364, -364, 128, 600}};