// code to build replacement list for gene labels based on pairwise BLAST output

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAXL 100             // maximum length (in characters) for a gene label
#define MAXLINE 1000         // maximum length for a line of BLAST output
#define MAXP 1000            // maximum number of partners in a gene label group
#define MAXN 10000           // maximum number of distinct gene labels

#define VERBOSE 0

// get index to a given string in a set of strings; return -1 if absent
// "labeldb" is a 1D array of characters; each of "nlabel" strings occupies MAXL characters in it
// we want n, where the string "label" we're searching for starts at n*MAXL (and ends at (n+1)*MAXL-1)
int getref(char *label, char *labeldb, int nlabel)
{
  int i;
  for(i = 0; i < nlabel; i++)
    {
      if(strcmp(label, &labeldb[MAXL*i]) == 0)
	return i;
    }
  return -1;
}

// placeholder to return BLAST quality from arguments
// this just performs a calculation and returns a value; it's given its own function to keep things tidy below
double score(int len1, int len2, int len12, float prop)
{
  double length;
  
  if((double)len12/len2 < (double)len12/len1)
    length = (double)len12/len2;
  else
    length = (double)len12/len1;

  return length*prop/100.;
}

int main(int argc, char *argv[])
{
  FILE *fp, *outputfp;
  int nlabel;
  int labelcount[MAXN];
  char *labeldb;
  char tstr[MAXLINE];
  char *partners;
  int npartners[MAXN];
  int addedpartners[MAXN];
  int ref1, ref2, pref1, pref2;
  int i, j, k;
  int jref, kref;
  int change;
  int starti;
  char label1[MAXL], label2[MAXL];
  float tmpf;
  double threshold = 0.25, sharedthreshold = 0.25;
  int len1, len2, len12;
  float prop;
  float thisscore;
  int *partnercount;
  int hitcount[MAXN];
  long int lcount;
  int bestcount, bestref;
  int done[MAXN];
  int iteration;
  
  // check that command-line arguments take the form we need
  if(argc < 6)
    {
      printf("Need label count file, BLAST output, output replacement file, score threshold, and intersection threshold\n");
      return 0;
    }

  // allocate memory for the arrays we'll be using
  // we'll use 1D arrays for everything
  labeldb = (char*)malloc(sizeof(char)*MAXL*MAXN);
  partners = (char*)malloc(sizeof(char)*MAXL*MAXN*MAXP);
  partnercount = (int*)malloc(sizeof(int)*MAXN*MAXN);

  threshold = atof(argv[4]);
  sharedthreshold = atof(argv[5]);

  /////
  // the first thing we do is read the occurrence counts for each label, produced by code upstream in the pipeline
  // this will help us decide which is the "most popular" label for a given group
  
  // read in labels and counts
  fp = fopen(argv[1], "r");
  if(fp == NULL)
    {
      printf("Couldn't open label count file %s\n", argv[1]);
      return 0;
    }

  // initialise gene label count
  nlabel = 0;

  // read and discard the csv header
  fgets(tstr, MAXL, fp);
  do{
    // read the next line; break if we've reached the end
    fgets(tstr, MAXL, fp);
    if(feof(fp)) break;
    // move through this line until we find a comma, storing the characters which make up the gene label
    for(i = 0; tstr[i] != ','; i++)
      labeldb[MAXL*nlabel + i] = tstr[i];
    // terminate gene label with a null character marking the end of the string
    labeldb[MAXL*nlabel + i] = '\0';
    // store the occurrence count which follows this point
    labelcount[nlabel] = atoi(&tstr[i+1]);
    // initialise partner count for this label, and increment the number found
    sprintf(&partners[MAXL*MAXP*nlabel], "%s", &labeldb[MAXL*nlabel]);
    npartners[nlabel] = 1;    
    nlabel++;
  }while(!feof(fp));
  fclose(fp);
  if(VERBOSE) printf("Read %i labels\n", nlabel);

  if(nlabel > MAXN)
    {
      printf("Need more memory allocated!\n");
      exit(1);
    }
  
  /////
  // next we go through the BLAST output, storing pairs where there's a convincing BLAST hit linking two different gene labels

  // go through BLAST output
  fp = fopen(argv[2], "r");
  if(fp == NULL)
    {
      printf("Couldn't open BLAST output file %s\n", argv[2]);
      return 0;
    }
  // initialise hit counter and line counter
  for(i = 0; i < MAXN; i++)
    hitcount[i] = 0;
  lcount = 0;
  do{
    // read the next line of BLAST output; break if we've reached the end
    fgets(tstr, MAXLINE, fp);
    if(feof(fp)) break;

    // update the user
    if(lcount % 1000000 == 0)
      printf("Line %li\n", lcount);
    lcount++;

    // (recall every line of the BLAST output is a pair of genes for which an overlap has been found; we want the gene labels and numerical features of the overlap)
    
    // get to the character following the first comma (that is, the start of the gene label, following the species)
    for(i = 0; !(tstr[i] == ','); i++);
    starti = i+1;
    // store the characters between here and the next comma as the first gene label
    for(i = starti; !(tstr[i] == ','); i++)
      label1[i-starti] = tolower(tstr[i]);
    // terminate gene label with a null character marking the end of the string
    label1[i-starti] = '\0';

    // get to the tab character that separates the two gene descriptions
    i += 2;
    for(; !(tstr[i] == '\t'); i++);
    // get to the character following the next comma (that is, the start of the gene label, following the second species)
    i++;
    for(; !(tstr[i] == ','); i++);
    starti = i+1;
    // store the characters between here and the next comma as the first gene label
    for(i = starti; !(tstr[i] == ','); i++)
      label2[i-starti] = tolower(tstr[i]);
    // terminate gene label with a null character marking the end of the string
    label2[i-starti] = '\0';

    // skip to the next tab, separating the numerical details of the overlap
    for(; tstr[i] != '\t'; i++);

    // if the two labels don't match, BLAST found a relationship between non-identical labels
    // otherwise, there's a relationship between two genes with the same label, which we ignore
    if(strcmp(label1, label2) && !(strcmp(label1, "anon\0") == 0 || strcmp(label2, "anon\0") == 0))
      {
	// get numerical info about the BLAST overlap -- read a value, skip to the next tab character, read the next value
	len1 = atoi(&tstr[i]); i++; for(;tstr[i] != '\t';i++); 
	len2 = atoi(&tstr[i]); i++; for(;tstr[i] != '\t';i++);
	len12 = atoi(&tstr[i]); i++; for(;tstr[i] != '\t';i++);
	prop = atof(&tstr[i]);

	// compute overlap score from these values
	thisscore = score(len1, len2, len12, prop);

	if(VERBOSE == 2) printf("Testing: %s %s %i %i %i %f = %f\n", label1, label2, len1, len2, len12, prop, thisscore);
	// if this is a "good" hit
	if(thisscore > threshold)
	  {
	    if(VERBOSE) printf("Seen a hit: %s %s %i %i %i %f = %f\n", label1, label2, len1, len2, len12, prop, thisscore);

	    // get indices for the two genes involved
	    // getref "looks up" the first argument in the second (which is an array concatenating N strings, where N is the third argument)
	    // it returns an index for the first argument if found, or -1 if not
	    ref1 = getref(label1, labeldb, nlabel);
	    ref2 = getref(label2, labeldb, nlabel);

	    if(ref1 == -1 || ref2 == -1)
	      {
		// this can be reached if the gene occurrence dataset doesn't contain one of the gene labels in the BLAST output
		// this means something's gone wrong with the upstream pipeline
		printf("Referencing error: couldn't find ");
		if(ref1 == -1) printf("%s ", label1);
  	        if(ref2 == -1) printf("%s ", label2);
		printf("\n");
		exit(1);
	      }

	    // check to see if gene2 is already stored as a partner of gene1
	    pref1 = getref(label2, &partners[MAXL*MAXP*ref1], npartners[ref1]);
	    if(pref1 == -1)
	      {
		// if not, add it to the array of partners
		// partners is a MAXN*MAXL*MAXP array. the string starting at partners[MAXL*MAXP*x + MAXL*y] is the yth partner of the label with reference x
		// npartners[x] is the current number of partners of the label with reference x
		// so this sprintf inserts "label2" as a new partner for the label with reference "ref1" (gene1)...
		printf("Adding %s as partner of %s (line %li)\n", &labeldb[MAXL*ref2], &labeldb[MAXL*ref1], lcount);
		sprintf(&partners[MAXL*MAXP*ref1 + MAXL*npartners[ref1]], "%s", label2);
		// this increments the number of *distinct* partners for gene1...
		npartners[ref1]++;
		// and this initialises the number of times this partner has been found for gene1
		partnercount[MAXN*ref1 + ref2] = 1;
		if(npartners[ref1] > MAXP)
		  {
		    printf("Need more partner memory!\n");
		  }
	      }
	    else
	      {
		// if gene2 is already a partner of gene1, increment the number of times it's been seen as such
  	        partnercount[MAXN*ref1 + ref2]++;
	      }
	    
	    // check to see if gene1 is already a partner of gene2
	    pref2 = getref(label1, &partners[MAXL*MAXP*ref2], npartners[ref2]);
	    if(pref2 == -1)
	      {
		printf("Adding %s as partner of %s\n", &labeldb[MAXL*ref1], &labeldb[MAXL*ref2]);
		// if not, add it; see comments immediately above
		sprintf(&partners[MAXL*MAXP*ref2 + MAXL*npartners[ref2]], "%s", label1);
		partnercount[MAXN*ref2 + ref1] = 1;
		npartners[ref2]++;
		if(npartners[ref2] > MAXP)
		  {
		    printf("Need more partner memory!\n");
		  }

	      }
	    else
	      partnercount[MAXN*ref2 + ref1]++;

	    hitcount[ref1]++;
	    hitcount[ref2]++;
	  }
      }
  }while(!feof(fp));
  fclose(fp);

  /////
  // now we prune pairs that look like "flukes" -- where the pair is only assigned based on a small subset (under "sharedthreshold") of appearances of one of the labels

  // loop through labels
  for(i = 0; i < nlabel; i++)
    {
      // loop through partners for this label
      for(j = 1; j < npartners[i]; j++)
	{
	  jref = getref(&partners[MAXL*MAXP*i + MAXL*j], labeldb, nlabel);
	  printf("Looking at %s -- %s (%i/%i) ", &labeldb[MAXL*i], &labeldb[MAXL*jref], partnercount[MAXN*jref + i], hitcount[jref]);
	  // if a small proportion of occurrences of this partner involve a link to label "i"
	  if(partnercount[MAXN*jref + i] < sharedthreshold*hitcount[jref])
	    {
	      printf(" -- chop\n");
	      // pop the partner from the list of "i"s partners
	      for(k = j; k < npartners[i]-1; k++)
		{
		  sprintf(&partners[MAXL*MAXP*i + MAXL*k], "%s", &partners[MAXL*MAXP*i + MAXL*(k+1)]);
		}
	      npartners[i]--;
	      j--;
	    }
	  else
	    {
	      printf(" -- keep\n");
	    }
	}
      // output list of "i"s partners
      if(npartners[i] > 1)
	{
	  for(j = 1; j < npartners[i]; j++)
	    printf("%s ", &partners[MAXL*MAXP*i + MAXL*j]);
	  printf("\n");
	}

    }
  
  /////
  // now we will fold together all linked labels, so that pairs are grouped together into sets
  // e.g if A-B and B-C, we want to form the A-C link too
  printf("Consolidating pairs...\n");
  iteration = 0;
  do{
    // here we're looping through each label i in our set
    // for each i, we loop through each partner j
    // for each j, we loop through THEIR partners k (ie next-nearest neighbours)
    // if k isn't a partner of i, we add it
    change = 0;
    for(i = 0; i < nlabel; i++)
      done[i] = 0;
    for(i = 0; i < nlabel; i++)
      {
	if(!done[i])
	  {
	    done[i] = 1;
	    
	    for(j = 1; j < npartners[i]; j++)
	      {
		jref = getref(&partners[MAXL*MAXP*i + MAXL*j], labeldb, nlabel);
		if(!done[jref])
		  {
		    printf("Looking at partner %s of %s\n", &partners[MAXL*MAXP*i + MAXL*j], &labeldb[MAXL*i]);

		    done[jref] = 1;
		    for(k = 1; k < npartners[jref]; k++)
		      {
			kref = getref(&partners[MAXL*MAXP*jref + MAXL*k], labeldb, nlabel);
			if(!done[kref])
			  {
			    done[kref] = 1;
			    printf("Looking at partner %s of partner %s of %s\n", &partners[MAXL*MAXP*jref + MAXL*k], &partners[MAXL*MAXP*i + MAXL*j], &labeldb[MAXL*i]); 
			    // check if k is a partner of i
			    pref1 = getref(&partners[MAXL*MAXP*jref + MAXL*k], &partners[MAXL*MAXP*i], npartners[i]);
			    if(pref1 == -1)
			      {
				// if not, add it and record that we've made a change
				sprintf(&partners[MAXL*MAXP*i + MAXL*npartners[i]], "%s", &partners[MAXL*MAXP*jref + MAXL*k]);
				partnercount[MAXN*i + npartners[i]] = 0;
				npartners[i]++;
				change = 1;
			      }
			  }
			/*		// check if i is a partner of k
					pref2 = getref(&labeldb[MAXL*i], &partners[MAXL*MAXP*kref], npartners[kref]);
					if(pref2 == -1)
					{
					// if not, add it and record that we've made a change
					sprintf(&partners[MAXL*MAXP*kref + MAXL*npartners[kref]], "%s", &labeldb[MAXL*i]);
					partnercount[MAXN*kref + npartners[kref]] = 0;
					npartners[kref]++;
					change = 1;
					}*/
		      }
		  }
	      }
      
	  }
      }
    iteration++;

  }while(change == 1); // keep looping until lists have converged
  
  
  outputfp = fopen(argv[3], "w");
  if(outputfp == NULL)
    {
      printf("Couldn't open output file %s\n", argv[3]);
      exit(1);
    }
  fprintf(outputfp, "OriginalLabel,ReplacementLabel\n");

  for(i = 0; i < MAXN; i++) 
    {
      done[i] = 0;
    }

  /////
  // finally, we ask which members of each set can reasonably be represented by the most common label in each set
  // this isn't every member, as some fluke overlaps may have occurred. we ask whether MANY (above "sharedthreshold") occurrences of a label are compatible with the most common label
  // this is done iteratively, so that e.g. if we have a group from A-B (reliable), B-C (reliable), B-D (fluke) with A most common we get:
  // first iteration (A-B); second iteration (A-B-C)
  // loop through found labels
  for(i = 0; i < nlabel; i++)
    {
      // if this label has partners and we haven't covered it yet
      if(npartners[i] > 1 && done[i] != 1)
	{
	  done[i] = 1;
	  printf("%s (%i) has %i partners: ", &labeldb[MAXL*i], labelcount[i], npartners[i]);
	  // loop through identified partners for each label, to find which is the most represented
	  bestcount = 0;
	  for(j = 0; j < npartners[i]; j++)
	    {
	      ref2 = getref(&partners[MAXL*MAXP*i + MAXL*j], labeldb, nlabel);
	      printf("%s ", &labeldb[MAXL*ref2]);
	      if(labelcount[ref2] > bestcount)
		{
		  // update most represented record
		  bestcount = labelcount[ref2];
		  bestref = ref2;
		}
	    }
	  printf(" -- most frequent is %s with %i\n", &labeldb[MAXL*bestref], bestcount);
	
	  // output set of replacements to file
	  for(j = 0; j < MAXL; j++)
	    {
	      addedpartners[j] = (j == bestref);
	    }
	  do{
	    change = 0;
	    // loop through "i"s partners
	    for(j = 0; j < npartners[i]; j++)
	      {
		jref = getref(&partners[MAXL*MAXP*i + MAXL*j], labeldb, nlabel);
		// if this partner hasn't been included in the replacement list under "i"s label yet
		if(!addedpartners[jref])
		  {
		    // loop through the list of partners that HAVE been included under "i"s label
		    for(k = 0; k < npartners[i] && !addedpartners[jref]; k++)
		      {
			kref = getref(&partners[MAXL*MAXP*i + MAXL*k], labeldb, nlabel);
			if(addedpartners[kref])
			  {
			    // if a good proportion of unincluded "j"s occurrences are paired with included "k"
			    if(partnercount[MAXN*jref + kref] > sharedthreshold*hitcount[jref])
			      {
				printf("  More than %.2f %s entries (%i/%i) are hits to %s, so adding to %s\n", sharedthreshold, &labeldb[MAXL*jref], partnercount[MAXN*jref + kref], hitcount[jref], &labeldb[MAXL*kref], &labeldb[MAXL*bestref]);
				fprintf(outputfp, "%s,%s\n", &labeldb[MAXL*jref], &labeldb[MAXL*bestref]);
				// include "j" in the label "i"
				addedpartners[jref] = 1;
				done[jref] = 1;
				change = 1;
			      }
			    else
			      {
				printf("  Fewer than %.2f %s entries (%i/%i) are hits to %s, so not yet adding to %s\n", sharedthreshold, &labeldb[MAXL*jref], partnercount[MAXN*jref + kref], hitcount[jref], &labeldb[MAXL*kref], &labeldb[MAXL*bestref]);
			      }
			  }
		      }
		  }
	      }
	  }while(change == 1);
	}
    }
  fclose(outputfp);

  return 0;
}
