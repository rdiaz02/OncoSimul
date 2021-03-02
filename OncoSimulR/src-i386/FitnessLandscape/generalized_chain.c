/*
	generalized chains
	only bi-allelic landscapes
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "landscape.h"
#include "genotypes.h"

void
int2index(struct landscape h, int i, int *k)
{
	int             l = h.ngenotypes;
	k[0] = i / l;
	k[1] = i % l;
}

int
index2int(struct landscape h, int *k)
{
	int             i;
	int             l = h.ngenotypes;
	i = (k[0]) * l + k[1];
	return i;
}

int
compare_genotypes(struct landscape h, int *genotype1, int *genotype2)
{
	//returns Hamming distance between genotypes
	int             i, j;
	int             l = h.nlocus;
	j = 0;
	for (i = 0; i < l; i++) {
		if (genotype1[i] != genotype2[i]) {
			j++;
		};
	};
	return j;
}

void
adjacency_matrix(struct landscape h, int *m)
{
	int             i, j, g, k[2];
	int             n = h.ngenotypes;
	int             l = h.nlocus;
	int             genotype1[l], genotype2[l];
	for (g = 0; g < n; g++) {
		int2genotype(h, g, genotype1);
		for (j = 0; j < n; j++) {
			k[0] = g;
			k[1] = j;
			i = index2int(h, k);
			if (h.fitness[j] <= h.fitness[g]) {
				m[i] = 0;
			} else {
				int2genotype(h, j, genotype2);
				if (compare_genotypes(h, genotype1, genotype2) == 1) {
					m[i] = 1;
				} else {
					m[i] = 0;
				};
			};
			printf("m= %u %u %u\n", g, j, m[i]);
		};
	};
}


void
initialize_flux(struct landscape h, double *flux)
{
	int             i;
	int             n = h.ngenotypes;
	for (i = 0; i < n * n; i++) {
		if (i % (n + 1) == 0) {
			flux[i] = 1.0;
		} else {
			flux[i] = -2.0;
		};
	};
}


double
get_flux(struct landscape h, int *m, int gi, int gf, double *flux)
{
	double          f;
	int             i, k[2], s;
	int             n = h.ngenotypes;
	k[1] = gf;
	k[0] = gi;
	if (flux[index2int(h, k)] < -1) {
		if (gi == gf) {
			f = 1.0;
		} else {
			if (h.fitness[gi] >= h.fitness[gf]) {
				f = 0.0;
			} else {
				f = 0;
				s = 0;
				k[0] = gi;
				for (i = 0; i < n; i++) {
					k[1] = i;
					if (m[index2int(h, k)] == 1) {
						//only use of the adjacency matrix m in this routine
							s++;
						f += get_flux(h, m, i, gf, flux);
					};
				};
				if (s > 0) {
					f = f / (double) s;
				} else {
					f = 0.0;
				};
			};
		};
		k[1] = gf;
		k[0] = gi;
		flux[index2int(h, k)] = f;
		printf("f= %u %u %f\n", gi, gf, f);
		return f;
	} else {
		return flux[index2int(h, k)];
	};
}

void
assign_fluxes(struct landscape h, int *m, double *flux)
{
	//finds the minima and calls get_flux from there
	int             i, j, s, k[2];
	int             n = h.ngenotypes;
	for (i = 0; i < n; i++) {
		s = 0;
		k[1] = i;
		for (j = 0; j < n; j++) {
			k[0] = j;
			s += m[index2int(h, k)];
		};
		if (s == 0) {
			printf("min %u\n", i);
			for (j = 0; j < n; j++) {
				k[0] = i;
				k[1] = j;
				flux[index2int(h, k)] = get_flux(h, m, i, j, flux);
			};
		};
	};
}


void
print_chain_pairs(struct landscape h, double *flux, double tol)
{
	int             i, j, k[2], c;
	int             l = h.nlocus;
	int             genotype[l];
	int             n = h.ngenotypes;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			k[0] = i;
			k[1] = j;
			if ((i != j) && (flux[index2int(h, k)] > 1 - tol)) {
				printf("%u\t%u\t", i, j);
				int2genotype(h, i, genotype);
				for (c = 0; c < l; c++)
					printf("%u", genotype[c]);
				printf("\t");
				int2genotype(h, j, genotype);
				for (c = 0; c < l; c++)
					printf("%u", genotype[c]);
				printf("\t%f", flux[index2int(h, k)]);
				printf("\n");
			};
		};
	};
}




//OTHER FUNCTIONS FOR NODE RELEVANCE ETC...


void
get_stats(struct landscape h, int g, int *m, int *np, int *nl, double *sigma, double *he, int fromall)
{
	int             j, s, u, k[2];
	int             n = h.ngenotypes;
	np[g] = 0;
	nl[g] = 0;
	sigma[g] = 0;
	s = 0;
	u = 0;
	k[0] = g;
	for (j = 0; j < n; j++) {
		k[1] = j;
		if (m[index2int(h, k)] == 1) {
			u++;
		};
	};
	if (u > 0) {
		he[g] = log((double) u);
	} else {
		he[g] = -1.0;
	};
	k[1] = g;
	for (j = 0; j < n; j++) {
		k[0] = j;
		if (m[index2int(h, k)] == 1) {
			s++;
			get_stats(h, j, m, np, nl, sigma, he, fromall);
			if (fromall == 0) {
				np[g] += np[j];
				nl[g] += (nl[j] + np[j]);
				sigma[g] += (sigma[j] + he[j] * np[j]);
			} else {
				np[g] += np[j] + 1;
				nl[g] += (nl[j] + np[j] + 1);
				sigma[g] += (sigma[j] + he[j] * (np[j] + 1));
			};
		};
	};
	if (s == 0) {
		printf("min %u\n", g);
		if (fromall == 0) {
			np[g] = 1;
		};
	};
}


void
initialize_stats(struct landscape h, int *m, int *np, int *nl, double *sigma, double *he, int fromall)
{

	int             i, j, s, k[2];
	int             n = h.ngenotypes;

	for (i = 0; i < n; i++) {

		s = 0;
		k[0] = i;

		for (j = 0; j < n; j++) {

			k[1] = j;
			s += m[index2int(h, k)];
		}

		if (s == 0) {

			printf("max %u\n", i);	/* COMPLETE !!!!!!!! */
			get_stats(h, i, m, np, nl, sigma, he, fromall);
		}
	}
}


void
get_node_constraints(struct landscape h, int *m, double *ha0, double *ha1)
{

	int             n = h.ngenotypes;
	int            *nl;
	int            *np;
	double         *sigma;
	double         *he;
	int             i;
	nl = malloc(sizeof(int) * n);
	np = malloc(sizeof(int) * n);
	sigma = malloc(sizeof(double) * n);
	he = malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) {
		nl[i] = 0;
		np[i] = 0;
		sigma[i] = 0;
		he[i] = 0;
	};
	initialize_stats(h, m, np, nl, sigma, he, 0);
	for (i = 0; i < n; i++) {
		if (nl[i] > 0) {
			ha0[i] = sigma[i] / (double) nl[i];
		} else {
			ha0[i] = 0;
		};
		nl[i] = 0;
		np[i] = 0;
		sigma[i] = 0;
		he[i] = 0;
	};
	initialize_stats(h, m, np, nl, sigma, he, 1);
	for (i = 0; i < n; i++) {
		if (nl[i] > 0) {
			ha1[i] = sigma[i] / (double) nl[i];
		} else {
			ha1[i] = 0;
		};
	};
	free(nl);
	free(np);
	free(sigma);
	free(he);
}
