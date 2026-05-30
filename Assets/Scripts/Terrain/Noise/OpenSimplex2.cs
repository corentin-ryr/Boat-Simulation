// Public-domain 2D OpenSimplex2 noise — direct port of KdotJPG's reference C# implementation
// (https://github.com/KdotJPG/OpenSimplex2). "Fast" variant: three lattice contributions per
// sample, output in roughly [-1, 1], isotropic (no axis-aligned banding like Perlin).
//
// Pure static functions, no shared mutable state — safe to call concurrently from worker
// threads. The static initializer fills the gradient table once; reads after that are
// allocation-free and lock-free.

namespace TerrainGrid.Noise
{
    public static class OpenSimplex2
    {
        // --- 2D constants (from KdotJPG's reference) ---

        const long PRIME_X = 0x5205402B9270C86FL;
        const long PRIME_Y = 0x598CD327003817B5L;
        const long HASH_MULTIPLIER = 0x53A3F72DEEC546F5L;

        const double SKEW_2D = 0.366025403784439;       // (sqrt(3) - 1) / 2
        const double UNSKEW_2D = -0.21132486540518713;  // (1/sqrt(3) - 1) / 2

        const int N_GRADS_2D_EXPONENT = 7;
        const int N_GRADS_2D = 1 << N_GRADS_2D_EXPONENT;
        const double NORMALIZER_2D = 0.05481866495625118;
        const float RSQUARED_2D = 2f / 3f;

        static readonly float[] GRADIENTS_2D;

        static OpenSimplex2()
        {
            // 24 base gradients evenly distributed around the unit circle (every 15°), normalized
            // so the noise output covers roughly [-1, 1] after the (R² - d²)⁴ falloff.
            float[] grad2 = {
                 0.38268343236509f,   0.923879532511287f,
                 0.923879532511287f,  0.38268343236509f,
                 0.923879532511287f, -0.38268343236509f,
                 0.38268343236509f,  -0.923879532511287f,
                -0.38268343236509f,  -0.923879532511287f,
                -0.923879532511287f, -0.38268343236509f,
                -0.923879532511287f,  0.38268343236509f,
                -0.38268343236509f,   0.923879532511287f,
                 0.130526192220052f,  0.99144486137381f,
                 0.608761429008721f,  0.793353340291235f,
                 0.793353340291235f,  0.608761429008721f,
                 0.99144486137381f,   0.130526192220052f,
                 0.99144486137381f,  -0.130526192220051f,
                 0.793353340291235f, -0.60876142900872f,
                 0.608761429008721f, -0.793353340291235f,
                 0.130526192220052f, -0.99144486137381f,
                -0.130526192220052f, -0.99144486137381f,
                -0.608761429008721f, -0.793353340291235f,
                -0.793353340291235f, -0.608761429008721f,
                -0.99144486137381f,  -0.130526192220052f,
                -0.99144486137381f,   0.130526192220051f,
                -0.793353340291235f,  0.608761429008721f,
                -0.608761429008721f,  0.793353340291235f,
                -0.130526192220052f,  0.99144486137381f,
            };

            for (int i = 0; i < grad2.Length; i++)
                grad2[i] = (float)(grad2[i] / NORMALIZER_2D);

            // Replicate to N_GRADS_2D entries (so the hash → gradient index is just a mask).
            GRADIENTS_2D = new float[N_GRADS_2D * 2];
            for (int i = 0, j = 0; i < GRADIENTS_2D.Length; i++, j++)
            {
                if (j == grad2.Length) j = 0;
                GRADIENTS_2D[i] = grad2[j];
            }
        }

        // 2D OpenSimplex2 noise. Output is approximately in [-1, 1].
        // (x, y) here use the standard 2D math convention; for terrain we pass (worldX, worldZ).
        public static float Noise2(long seed, double x, double y)
        {
            // Get points for A2* lattice (skewed simplex grid).
            double s = SKEW_2D * (x + y);
            double xs = x + s, ys = y + s;
            return Noise2_UnskewedBase(seed, xs, ys);
        }

        static float Noise2_UnskewedBase(long seed, double xs, double ys)
        {
            // Skewed simplex coords; find the containing simplex cell.
            int xsb = FastFloor(xs), ysb = FastFloor(ys);
            float xi = (float)(xs - xsb), yi = (float)(ys - ysb);
            long xsbp = xsb * PRIME_X, ysbp = ysb * PRIME_Y;
            float t = (xi + yi) * (float)UNSKEW_2D;
            float dx0 = xi + t, dy0 = yi + t;

            float value = 0f;

            // First vertex (origin of the containing cell).
            float a0 = RSQUARED_2D - dx0 * dx0 - dy0 * dy0;
            if (a0 > 0)
                value = (a0 * a0) * (a0 * a0) * Grad(seed, xsbp, ysbp, dx0, dy0);

            // Diagonally-opposite vertex.
            float a1 = (float)(2 * (1 + 2 * UNSKEW_2D) * (1 / UNSKEW_2D + 2)) * t
                     + ((float)(-2 * (1 + 2 * UNSKEW_2D) * (1 + 2 * UNSKEW_2D)) + a0);
            if (a1 > 0)
            {
                float dx1 = dx0 - (float)(1 + 2 * UNSKEW_2D);
                float dy1 = dy0 - (float)(1 + 2 * UNSKEW_2D);
                value += (a1 * a1) * (a1 * a1) * Grad(seed, xsbp + PRIME_X, ysbp + PRIME_Y, dx1, dy1);
            }

            // Third vertex — depends on which of the two triangles in the cell we're in.
            if (dy0 > dx0)
            {
                float dx2 = dx0 - (float)UNSKEW_2D;
                float dy2 = dy0 - (float)(UNSKEW_2D + 1);
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0)
                    value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp, ysbp + PRIME_Y, dx2, dy2);
            }
            else
            {
                float dx2 = dx0 - (float)(UNSKEW_2D + 1);
                float dy2 = dy0 - (float)UNSKEW_2D;
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0)
                    value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp + PRIME_X, ysbp, dx2, dy2);
            }

            return value;
        }

        static float Grad(long seed, long xsvp, long ysvp, float dx, float dy)
        {
            long hash = seed ^ xsvp ^ ysvp;
            hash *= HASH_MULTIPLIER;
            hash ^= hash >> (64 - N_GRADS_2D_EXPONENT - 1);
            int gi = (int)hash & ((N_GRADS_2D - 1) << 1);
            return GRADIENTS_2D[gi | 0] * dx + GRADIENTS_2D[gi | 1] * dy;
        }

        static int FastFloor(double x)
        {
            int xi = (int)x;
            return x < xi ? xi - 1 : xi;
        }
    }
}
