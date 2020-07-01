#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mem2/bwa.h"

static int
usage (void)
{
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:  fast_aligner index [-p prefix] <in.fasta>\n");
	fprintf (stderr, "\n");

	return 1;
}

int
pipe_db_build_main (int argc, char * argv[])
{
	if (argc < 2)
		return usage ();

	char * prefix;
	int c;
  int ret;

	prefix = NULL;
	while ((c = getopt(argc,argv,"p:")) != -1) {
		if (c == 'p')
			prefix = optarg;
		else
			return 1;
	}

	if (optind > argc-1)
		return usage ();

	if (prefix == NULL)
		prefix = argv[optind];

  // build FM-index for mapping
	bwa_idx_build (argv[optind], prefix);

  return 0;
}
