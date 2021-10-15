#include "msa.h"
#include "io.h"
#include "function.h"
#include <linux/limits.h>

#define REPORTCOSTS 1
#if REPORTCOSTS
#include <time.h>
#endif
#define VERSION "0.3.0"
#define SHOWVERSION reporterr( "%s (%s, %d-bit) Version " VERSION "\n\n", "MSA align", (seq_type == 1) ? "nuc" : ((seq_type == 0) ? "unknown" : "aa"), sizeof(int *) * 8 )
// #define FILESAVE

static int alignmode, calcsp;

void print_help_message()
{
    reporterr("WMSA alignment Version %s help: \n", VERSION);
    reporterr("== Common ==\n");
    reporterr("-i: input file name\n");
    reporterr("-o: output file name\n");
    reporterr("-D, -P: The sequence is DNA/RNA(-D) or Protein(-P)\n");
    reporterr("-T: use threads in cd-hit and staralign\n");
    reporterr("-p: Calcuate SP Scores after alignment\n");
    reporterr("-d: print debug info on profilealign, staralign and cd-hit*\n");
    reporterr("== MAFFT common arguments ==\n");
    reporterr("-V, -f, -S: ppenalty_dist, penalty, nmax_shift\n");
    reporterr("-z, -w: fftthreshold, fftWinsize\n");
    reporterr("-b: Band in alignment\n");
    reporterr("-B: BLOSUM??\n");
    reporterr("=== profilealign arguments ===\n");
    reporterr("-A, -F: Use Simply lign(-A) or FFT align(-F) on profile alignment\n");
    reporterr("== CD-HIT arugments ==\n");
    reporterr("-c: threshold on cd-hit to cluster\n");
    reporterr("-M: max memory for CD-HIT\n");
}

void arguments(int argc, char *argv[])
{
    if(argc == 1) 
    {
        print_help_message();
        exit(0);
    }
    int c;
    cdhitsim = NOTKNOWNDOUBLE;
    ppenalty = NOTKNOWNINT;
    ppenalty_dist = NOTKNOWNINT;
    fftthreshold = NOTKNOWNINT;
    fftWinSize = NOTKNOWNINT;
    alignband = 20;
    threads = 1;
    nmax_shift = 1;
    seq_type = 0;
    alignmode = 1;
    calcsp = 0;
    maxmemory = 1024;
    printdebug = 0;
    BLOSUM = 62;
    while(--argc > 0 && (*++ argv)[0] == '-' )
	{
        while ( (c = *++ argv[0]) )
		{
            switch( c )
            {
				// I/O
				case 'i':
					inputfile = *++ argv;
					-- argc;
					goto nextoption;
                case 'o':
                    outputfile = *++ argv;
                    -- argc;
                    goto nextoption;
                // Sequence type
                case 'P':
                    seq_type = 2;
                    goto nextoption;
                case 'D':
                    seq_type = 1;
                    goto nextoption;
				// MAFFT arguments
				case 'V':
					ppenalty_dist = myatof( *++ argv );
					-- argc;
					goto nextoption;
				case 'f':
					ppenalty = myatof( *++ argv );
					-- argc;
					goto nextoption;
				case 'z':
					fftthreshold = myatoi( *++ argv );
					-- argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++ argv );
					-- argc;
					goto nextoption;
                case 'b':
                    alignband = myatoi(*++ argv);
                    -- argc;
                    goto nextoption;
                case 'T':
                    threads = myatoi(*++ argv);
                    -- argc;
                    goto nextoption;
                case 'S':
                    nmax_shift = myatoi(*++ argv);
                    -- argc;
                    goto nextoption;
                case 'B':
                    BLOSUM = myatoi(*++ argv);
                    -- argc;
                    goto nextoption;
                case 'A':
                    alignmode = 0;
                    break;
                case 'F':
                    alignmode = 1;
                    break;
                case 'd':
                    printdebug = 1;
                    break;
				// CD-HIT arguments
				case 'c':
					cdhitsim = myatof( *++ argv );
                    if(cdhitsim >= 1.0 || cdhitsim <= 0.0) 
                    {
                        reporterr("Warning: cd-hit-sim is out of (0, 1), will ignore. \n");
                        cdhitsim = NOTKNOWNDOUBLE;
                    }
					-- argc;
					goto nextoption;
                case 'M':
                    maxmemory = myatoi(*++ argv);
                    -- argc;
                    goto nextoption;
                // common
                case 'p':
                    calcsp = 1;
                    break;
                case 'H':
                case '?':
                    print_help_message();
                    exit(0);
                default:
                    reporterr( "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc != 0 ) 
    {
        reporterr( "Please type wmsa -H or -? to get help.\n" );
        exit( 1 );
    }
}

int main(int argc, char **argv) 
{
    arguments(argc, argv);
    FILE *input = fopen(inputfile, "rb"), *cmdd, *order;
    seq_type = getnumlenandtype(input, seq_type);
    reporterr("This file has %d sequences.\n", seq_num);
    
    /* Part -1: strip gap */
    char **aseq, **bseq, **strip_name;
    int *nlen;
    int i, j, k;
    aseq = AllocateCharMtx(seq_num, seq_maxlen + 10);
    bseq = AllocateCharMtx(seq_num, seq_maxlen + 10);
    strip_name = AllocateCharMtx(seq_num, B + 1);
    nlen = AllocateIntVec(seq_num + 1);
    readData(input, strip_name, nlen, aseq, seq_num);
    fclose(input);
    if(seq_num < 2)
    {
        reporterr("Warning: only %d sequence(s) found.\n", seq_num);
        return 0;
    }
    for(i = 0; i < seq_num; ++ i) gapfilter_oneseq(bseq[i], aseq[i]);
    cmdd = fopen(outputfile, "w");
    writeData(cmdd, seq_num, strip_name, nlen, bseq);
    FreeCharMtx(aseq);
    FreeCharMtx(bseq);
    FreeCharMtx(strip_name);
    FreeIntVec(nlen);
    fclose(cmdd);
    
#if REPORTCOSTS
		time_t starttime;
		starttime = time(NULL);
#endif
    /* Part 0: process program path && make temp folder */
    // get the program folder
    char *programfolder = AllocateCharVec(PATH_MAX + 100);
    programfolder = get_exe_path(programfolder, PATH_MAX + 100);
    if(printdebug) puts(programfolder);
    // now the value programfolder is the right value, and the last place of the folder name is not '/'.
    if(maketmpfolder())
    {
        reporterr("Warning: cannot make folder %s\n", tmpdir);
    }
    /* Part 1: cd-hit process input */
    reporterr("Clustering...");
    clustercommand(outputfile, programfolder); 
    reporterr("done. \n");
    /* 
       Saved file to tmpdir/tmp;
       the clusters is tmpdir/tmp_%d.clstr
    */

    /* Part 2: for each cluster, run staralign to align */
    int cluster_seq = calcuate_cluster_sequences();
    reporterr("CD-HIT spilted this file to %d cluster(s).\n", cluster_seq);
    
    /* Part 2.1: read the center file */
    int *center_len;
    char **center_name = NULL, **center_seq = NULL;

    /* Part 2.2: staralign on the code && write center sequence order file */
    reporterr("BLASTalign...\n");
    // first check the existance of this program
    checkBLASTalign(programfolder, "mafft/staralign");
    for(i = 0; i < cluster_seq; ++ i)
    {
        sprintf(cmdstr2, "%stmp_%d.clstr", tmpdir, i);
        sprintf(cmdstr3, "%stmp_%d.center", tmpdir, i);
        BLASTaligncommand(cmdstr3, cmdstr2, programfolder, "mafft/staralign");
        reporterr("\rSTEP %d / %d", i + 1, cluster_seq);
    }
    /* Part 2.3: Write center sequence order file */
    sprintf(cmdstr, "%scluster_order", tmpdir);
    cmdd = fopen(cmdstr, "w");
    for(i = 0; i < cluster_seq; ++ i)
    {
        fprintf(cmdd, "%stmp_%d.clstr.res\n", tmpdir, i);
    }
    fclose(cmdd);
    reporterr("\ndone. \n");
    
    /* Part 3: profile-profile align */
    reporterr("profile merging... ");
    checkprofilealign(programfolder, "mafft/profilealign");
    sprintf(cmdstr3, "%scluster_order", tmpdir);
    sprintf(cmdstr2, "%stmp", tmpdir);
    profilealigncommand(cmdstr3, cmdstr2, alignmode, programfolder, "mafft/profilealign");
    reporterr("done. \n");

#if REPORTCOSTS
    reporterr( "\nmsa align, real = %f min\n", (float)(time(NULL) - starttime)/60.0 );
#endif

    reporterr("\n\nDone.\n");
    SHOWVERSION;
    if(calcsp)
    {
        /* SP Scores calcuation */
        if(seq_type == 2) reporterr("Info: Using BLOSUM62 to calcuate the SP Scores\n");
        else if(seq_type == 1) reporterr("Info: Using trans/trans model to calcuate the SP Scores\n");
        input = fopen(outputfile, "rb");
        getnumlenandtype(input, seq_type);
        reporterr("Info: The final file has %d sequences.\n", seq_num);
        center_name = AllocateCharMtx(seq_num, seq_maxlen + 10);
        center_len = AllocateIntVec(seq_num);
        center_seq = AllocateCharMtx(seq_num, seq_maxlen + 10);
        readData(input, center_name, center_len, center_seq, seq_num);
        fclose(input);
        double sp = 0.0;
        if(ppenalty == NOTKNOWNINT) ppenalty = (int)(-1530.0 * 600 / 1000 + 0.5); // by MAFFT
        
        for(i = 0; i < seq_num; ++ i)
            for(j = i + 1; j < seq_num; ++ j)
            {
                sp += naivepairscore11(center_seq[i], center_seq[j], ppenalty, seq_type) / 600.0;
            }
        reporterr("SP Score is %f\n", sp);
        FreeCharMtx(center_name);
        FreeCharMtx(center_seq);
        FreeIntVec(center_len);
    }
#ifndef FILESAVE
    /* Part 4: remove useless files */
    cleantmpfile();
#endif
    free(programfolder);
    return 0;
}