import java.util.Arrays;

public class SchwarzschildEinstein {

    static final int DIM = 4; // (t, r, theta, phi)

    public static void main(String[] args) {
        double M = 1.0;                 // masse (unités géométrisées)
        double r = 10.0;                // point d'évaluation
        double theta = Math.PI / 3.0;   // point d'évaluation
        double t = 0.0;
        double phi = 0.0;

        double[] x = new double[]{t, r, theta, phi};

        // Pas de dérivation (à ajuster si besoin)
        double hr = 1e-5 * Math.max(1.0, r);
        double hth = 1e-6;
        double[] h = new double[]{1e-3, hr, hth, 1e-3}; // t, r, theta, phi

        double[][] G = einsteinTensorAt(x, h, M);

        System.out.println("G_{mu nu} at (r,theta)=(" + r + "," + theta + "), M=" + M);
        printMatrix(G);

        // En Schwarzschild vide, on s'attend à ~0 (erreurs numériques petites)
        System.out.println("\nMax |G_{mu nu}| = " + maxAbs(G));
    }

    // ---------------------------
    // Schwarzschild metric g_{mu nu}
    // ---------------------------
    static double[][] metric(double[] x, double M) {
        double r = x[1];
        double th = x[2];

        double f = 1.0 - 2.0 * M / r;

        double[][] g = new double[DIM][DIM];
        // diag(-f, 1/f, r^2, r^2 sin^2 th)
        g[0][0] = -f;
        g[1][1] = 1.0 / f;
        g[2][2] = r * r;
        double s = Math.sin(th);
        g[3][3] = r * r * s * s;
        return g;
    }

    // ---------------------------
    // Einstein tensor G_{mu nu}
    // ---------------------------
    static double[][] einsteinTensorAt(double[] x, double[] h, double M) {
        double[][] g = metric(x, M);
        double[][] gInv = invert4x4(g);

        // dg[a][m][n] = ∂_a g_{mn}
        double[][][] dg = metricPartialsCentral(x, h, M);

        // Gamma[rho][mu][nu] = Γ^rho_{mu nu}
        double[][][] Gamma = christoffel(gInv, dg);

        // dGamma[a][rho][mu][nu] = ∂_a Γ^rho_{mu nu} via différences finies
        double[][][][] dGamma = christoffelPartialsCentral(x, h, M);

        // Riemann: R[rho][sigma][mu][nu] = R^rho_{ sigma mu nu}
        double[][][][] Riem = riemann(Gamma, dGamma);

        // Ricci: Ric[mu][nu] = R^rho_{ mu rho nu}
        double[][] Ric = ricci(Riem);

        // scalaire: R = g^{mu nu} Ric_{mu nu}
        double Rscalar = 0.0;
        for (int mu = 0; mu < DIM; mu++) {
            for (int nu = 0; nu < DIM; nu++) {
                Rscalar += gInv[mu][nu] * Ric[mu][nu];
            }
        }

        // Einstein: G_{mu nu} = Ric_{mu nu} - 1/2 g_{mu nu} R
        double[][] G = new double[DIM][DIM];
        for (int mu = 0; mu < DIM; mu++) {
            for (int nu = 0; nu < DIM; nu++) {
                G[mu][nu] = Ric[mu][nu] - 0.5 * g[mu][nu] * Rscalar;
            }
        }
        return G;
    }

    // ---------------------------
    // Metric partial derivatives (central differences)
    // ---------------------------
    static double[][][] metricPartialsCentral(double[] x, double[] h, double M) {
        double[][][] dg = new double[DIM][DIM][DIM];
        for (int a = 0; a < DIM; a++) {
            double ha = h[a];
            if (ha == 0) continue;

            double[] xp = x.clone();
            double[] xm = x.clone();
            xp[a] += ha;
            xm[a] -= ha;

            double[][] gp = metric(xp, M);
            double[][] gm = metric(xm, M);

            for (int m = 0; m < DIM; m++) {
                for (int n = 0; n < DIM; n++) {
                    dg[a][m][n] = (gp[m][n] - gm[m][n]) / (2.0 * ha);
                }
            }
        }
        return dg;
    }

    // ---------------------------
    // Christoffel Γ^rho_{mu nu} from g^{-1} and ∂g
    // Γ^rho_{mu nu} = 1/2 g^{rho sigma} (∂_mu g_{sigma nu} + ∂_nu g_{sigma mu} - ∂_sigma g_{mu nu})
    // ---------------------------
    static double[][][] christoffel(double[][] gInv, double[][][] dg) {
        double[][][] Gamma = new double[DIM][DIM][DIM];
        for (int rho = 0; rho < DIM; rho++) {
            for (int mu = 0; mu < DIM; mu++) {
                for (int nu = 0; nu < DIM; nu++) {
                    double sum = 0.0;
                    for (int sig = 0; sig < DIM; sig++) {
                        sum += gInv[rho][sig] * (dg[mu][sig][nu] + dg[nu][sig][mu] - dg[sig][mu][nu]);
                    }
                    Gamma[rho][mu][nu] = 0.5 * sum;
                }
            }
        }
        return Gamma;
    }

    // ---------------------------
    // ∂_a Γ^rho_{mu nu} (central differences, recompute Γ)
    // ---------------------------
    static double[][][][] christoffelPartialsCentral(double[] x, double[] h, double M) {
        double[][][][] dGamma = new double[DIM][DIM][DIM][DIM];

        for (int a = 0; a < DIM; a++) {
            double ha = h[a];
            if (ha == 0) continue;

            double[] xp = x.clone();
            double[] xm = x.clone();
            xp[a] += ha;
            xm[a] -= ha;

            double[][] gp = metric(xp, M);
            double[][] gm = metric(xm, M);

            double[][] gInvP = invert4x4(gp);
            double[][] gInvM = invert4x4(gm);

            double[][][] dgP = metricPartialsCentral(xp, h, M);
            double[][][] dgM = metricPartialsCentral(xm, h, M);

            double[][][] GammaP = christoffel(gInvP, dgP);
            double[][][] GammaM = christoffel(gInvM, dgM);

            for (int rho = 0; rho < DIM; rho++) {
                for (int mu = 0; mu < DIM; mu++) {
                    for (int nu = 0; nu < DIM; nu++) {
                        dGamma[a][rho][mu][nu] = (GammaP[rho][mu][nu] - GammaM[rho][mu][nu]) / (2.0 * ha);
                    }
                }
            }
        }
        return dGamma;
    }

    // ---------------------------
    // Riemann: R^rho_{ sigma mu nu}
    // = ∂_mu Γ^rho_{nu sigma} - ∂_nu Γ^rho_{mu sigma}
    //   + Γ^rho_{mu λ} Γ^λ_{nu σ} - Γ^rho_{nu λ} Γ^λ_{mu σ}
    // ---------------------------
    static double[][][][] riemann(double[][][] Gamma, double[][][][] dGamma) {
        double[][][][] R = new double[DIM][DIM][DIM][DIM];
        for (int rho = 0; rho < DIM; rho++) {
            for (int sig = 0; sig < DIM; sig++) {
                for (int mu = 0; mu < DIM; mu++) {
                    for (int nu = 0; nu < DIM; nu++) {
                        double val = dGamma[mu][rho][nu][sig] - dGamma[nu][rho][mu][sig];
                        for (int lam = 0; lam < DIM; lam++) {
                            val += Gamma[rho][mu][lam] * Gamma[lam][nu][sig];
                            val -= Gamma[rho][nu][lam] * Gamma[lam][mu][sig];
                        }
                        R[rho][sig][mu][nu] = val;
                    }
                }
            }
        }
        return R;
    }

    // ---------------------------
    // Ricci: R_{mu nu} = R^rho_{ mu rho nu}
    // ---------------------------
    static double[][] ricci(double[][][][] Riem) {
        double[][] Ric = new double[DIM][DIM];
        for (int mu = 0; mu < DIM; mu++) {
            for (int nu = 0; nu < DIM; nu++) {
                double sum = 0.0;
                for (int rho = 0; rho < DIM; rho++) {
                    sum += Riem[rho][mu][rho][nu];
                }
                Ric[mu][nu] = sum;
            }
        }
        return Ric;
    }

    // ---------------------------
    // 4x4 matrix inversion (Gauss-Jordan)
    // ---------------------------
    static double[][] invert4x4(double[][] A) {
        int n = 4;
        double[][] aug = new double[n][2 * n];

        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, aug[i], 0, n);
            aug[i][n + i] = 1.0;
        }

        for (int col = 0; col < n; col++) {
            // pivot
            int pivot = col;
            double max = Math.abs(aug[col][col]);
            for (int row = col + 1; row < n; row++) {
                double v = Math.abs(aug[row][col]);
                if (v > max) {
                    max = v;
                    pivot = row;
                }
            }
            if (max < 1e-18) {
                throw new IllegalArgumentException("Matrix seems singular (pivot too small).");
            }
            if (pivot != col) {
                double[] tmp = aug[col];
                aug[col] = aug[pivot];
                aug[pivot] = tmp;
            }

            // normalize pivot row
            double p = aug[col][col];
            for (int j = 0; j < 2 * n; j++) aug[col][j] /= p;

            // eliminate other rows
            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double factor = aug[row][col];
                for (int j = 0; j < 2 * n; j++) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        double[][] inv = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(aug[i], n, inv[i], 0, n);
        }
        return inv;
    }

    // ---------------------------
    // Utilities
    // ---------------------------
    static void printMatrix(double[][] A) {
        for (double[] row : A) {
            System.out.println(Arrays.toString(row));
        }
    }

    static double maxAbs(double[][] A) {
        double m = 0.0;
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                m = Math.max(m, Math.abs(A[i][j]));
            }
        }
        return m;
    }
}
