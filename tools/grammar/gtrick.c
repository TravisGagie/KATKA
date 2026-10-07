// Size of the SPIRE'20 encoding (Gagie et al., "Practical random access to SLP-compressed texts") of a bigrepair -i grammar:
// non-terminals grouped by expansion length; each rule stored as (left child's length, left child's offset in its group,
// right child's offset in its group); terminals: length 1, offset = symbol.  Start rule: lengths in unary (Elias-Fano) + offsets.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static int lg(uint64_t x) { int b = 0; while ((1ULL << b) < x) b++; return b; }   // ceil(log2 x)
int main(int argc, char **argv) {
  char fn[4096]; sprintf(fn, "%s.R", argv[1]); FILE *f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  int alpha; if (fread(&alpha, 4, 1, f) != 1) return 1; long R = (sz - 4) / 8; int32_t *rl = malloc(8 * R); if (fread(rl, 8, R, f) != (size_t)R) return 1; fclose(f);
  long S = alpha + R; uint32_t *len = malloc(4 * S); for (int i = 0; i < alpha; i++) len[i] = 1;
  uint32_t maxl = 1; for (long i = 0; i < R; i++) { len[alpha+i] = len[rl[2*i]] + len[rl[2*i+1]]; if (len[alpha+i] > maxl) maxl = len[alpha+i]; }
  uint64_t *gsz = calloc(maxl + 1, 8); for (long i = 0; i < R; i++) gsz[len[alpha+i]]++;
  gsz[1] = alpha;   // terminals form the length-1 group
  long d = 0; for (uint32_t l = 2; l <= maxl; l++) if (gsz[l]) d++;
  sprintf(fn, "%s.C", argv[1]); f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); long C = ftell(f) / 4; fseek(f, 0, SEEK_SET);
  int32_t *c = malloc(4 * C); if (fread(c, 4, C, f) != (size_t)C) return 1; fclose(f);
  // fixed-width per group: offsets of a child in a group of size s cost ceil(log2 s) bits; left length costs ceil(log2(d+1)) bits
  double bl = lg(d + 2), boff = 0, bc = 0; uint64_t N = 0;
  for (long i = 0; i < R; i++) { boff += lg(gsz[len[rl[2*i]]] + 1) + lg(gsz[len[rl[2*i+1]]] + 1); }
  for (long i = 0; i < C; i++) { bc += lg(gsz[len[c[i]]] + 1); N += len[c[i]]; }
  double ef = C * (2 + log2((double)N / C));   // Elias-Fano for the start rule's lengths
  double mphf = 1.56 * d;
  double bl2 = 0; for (long i = 0; i < R; i++) bl2 += lg(len[alpha+i] - 1); printf("  left length within its group (ceil log2(l-1) bits): %.1f bits per rule\n", bl2 / R);
  printf("  variant total %.3f GB\n", (bl2 + boff + bc + ef + mphf) / 8e9);
  double tot = R * bl + boff + bc + ef + mphf;
  printf("%s: rules %ld, distinct expansion lengths d %ld, |C| %ld, n %lu\n", argv[1], R, d, C, N);
  printf("  per rule: left length %.0f bits + offsets %.1f bits = %.1f bits\n", bl, boff / R, bl + boff / R);
  printf("  start rule: offsets %.1f bits/symbol, Elias-Fano lengths %.3f MB\n", bc / C, ef / 8e6);
  printf("  SPIRE'20 encoding total %.3f GB (lengths included)\n", tot / 8e9);
  printf("  naive: 2 x %d bits per rule + start rule = %.3f GB, plus %d-bit lengths = %.3f GB\n", lg(S), (2.0 * R * lg(S) + C * lg(S)) / 8e9, lg(N), (2.0 * R * lg(S) + C * lg(S) + R * (double)lg(N)) / 8e9);
  return 0;
}
