// Cutting a bigrepair -i grammar at expansion length w: symbols with expansion < w whose parent (or the start rule) has
// expansion >= w are stored explicitly (their expansions, in a dictionary); rules with expansion >= w stay as rules.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int lg(uint64_t x) { int b = 0; while ((1ULL << b) < x) b++; return b; }
int main(int argc, char **argv) {
  char fn[4096]; sprintf(fn, "%s.R", argv[1]); FILE *f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  int alpha; if (fread(&alpha, 4, 1, f) != 1) return 1; long R = (sz - 4) / 8; int32_t *rl = malloc(8 * R); if (fread(rl, 8, R, f) != (size_t)R) return 1; fclose(f);
  long S = alpha + R; uint32_t *len = malloc(4 * S); for (int i = 0; i < alpha; i++) len[i] = 1;
  for (long i = 0; i < R; i++) len[alpha+i] = len[rl[2*i]] + len[rl[2*i+1]];
  sprintf(fn, "%s.C", argv[1]); f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); long C = ftell(f) / 4; fseek(f, 0, SEEK_SET);
  int32_t *c = malloc(4 * C); if (fread(c, 4, C, f) != (size_t)C) return 1; fclose(f);
  uint8_t *cut = malloc(S); int ws[] = {4, 8, 16, 32, 64, 128, 256};
  printf("w  upper_rules  cut_symbols  dict_tokens  dict_GB(14b)  upper_GB  total_GB(+EF starts 0.228)\n");
  for (int t = 0; t < 7; t++) {
    uint32_t w = ws[t]; long up = 0; for (long s = 0; s < S; s++) cut[s] = 0;
    for (long i = 0; i < R; i++) if (len[alpha+i] >= w) { up++; for (int k = 0; k < 2; k++) { int ch = rl[2*i+k]; if (len[ch] < w) cut[ch] = 1; } }
    for (long i = 0; i < C; i++) if (len[c[i]] < w) cut[c[i]] = 1;
    long nc = 0; uint64_t dt = 0; for (long s = 0; s < S; s++) if (cut[s]) { nc++; dt += len[s]; }
    double idb = lg(up + nc + 1);   // a symbol of the upper part or the start rule names an upper rule or a cut symbol
    double upper = (2.0 * up + C) * idb / 8e9, dict = dt * 14.0 / 8e9 + nc * (lg(dt) ) / 8e9;   // + offsets into the dictionary
    printf("%-3u %11ld %12ld %12lu %13.3f %9.3f %9.3f\n", w, up, nc, dt, dict, upper, dict + upper + 0.228);
  }
  return 0;
}
