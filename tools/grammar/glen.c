// Expansion lengths (in tokens = runs) of a RePair grammar from bigrepair -i (.R: alpha, then rule pairs; .C: sequence).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int cmp(const void *a, const void *b) { uint64_t x = *(uint64_t*)a, y = *(uint64_t*)b; return x < y ? -1 : x > y; }
int main(int argc, char **argv) {
  char fn[4096]; sprintf(fn, "%s.R", argv[1]); FILE *f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  int alpha; fread(&alpha, 4, 1, f); long R = (sz - 4) / 8; int32_t *rl = malloc(8 * R); fread(rl, 8, R, f); fclose(f);
  uint32_t *len = malloc(4 * (alpha + R)); for (int i = 0; i < alpha; i++) len[i] = 1;
  int maxd = 0; uint8_t *dep = malloc(alpha + R); for (int i = 0; i < alpha; i++) dep[i] = 0;
  for (long i = 0; i < R; i++) { int a = rl[2*i], b = rl[2*i+1]; len[alpha+i] = len[a] + len[b]; int d = 1 + (dep[a] > dep[b] ? dep[a] : dep[b]); dep[alpha+i] = d > 255 ? 255 : d; if (d > maxd) maxd = d; }
  sprintf(fn, "%s.C", argv[1]); f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); long C = ftell(f) / 4; fseek(f, 0, SEEK_SET);
  int32_t *c = malloc(4 * C); fread(c, 4, C, f); fclose(f);
  uint64_t tot = 0, *cl = malloc(8 * C); int cmaxd = 0; long cdsum = 0;
  for (long i = 0; i < C; i++) { cl[i] = len[c[i]]; tot += cl[i]; cdsum += dep[c[i]]; if (dep[c[i]] > cmaxd) cmaxd = dep[c[i]]; }
  printf("%s: alpha %d, rules %ld, |C| %ld, expanded length %lu, max rule depth %d (capped 255), mean depth of C symbols %.1f\n", argv[1], alpha, R, C, tot, maxd, (double)cdsum / C);
  qsort(cl, C, 8, cmp);
  printf("expansion of C symbols: mean %.1f, median %lu, p90 %lu, p99 %lu, p99.9 %lu, max %lu\n", (double)tot / C, cl[C/2], cl[(long)(.9*C)], cl[(long)(.99*C)], cl[(long)(.999*C)], cl[C-1]);
  long thr[] = {64, 256, 1024, 4096, 16384};
  for (int t = 0; t < 5; t++) { long k = 0; double bits = 0; for (long i = 0; i < R; i++) if (len[alpha+i] >= thr[t]) { k++; } 
    printf("rules with expansion >= %ld: %ld (%.2f%%)\n", thr[t], k, 100.0 * k / R); }
  // bits for storing every rule's length with an Elias-gamma-like variable-length code
  double gb = 0; for (long i = 0; i < R; i++) { uint32_t x = len[alpha+i]; int lg = 31 - __builtin_clz(x); gb += 2 * lg + 1; }
  printf("all rule lengths: %.1f bits each with Elias gamma (%.3f GB), vs 29 bits fixed (%.3f GB)\n", gb / R, gb / 8e9, R * 29 / 8e9);
  return 0;
}
