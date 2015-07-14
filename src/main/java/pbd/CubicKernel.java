package pbd;

import com.badlogic.gdx.math.MathUtils;
import com.badlogic.gdx.math.Vector3;

public class CubicKernel {
	private static final float K_FACTOR = 8.0f;
	private static final float L_FACTOR = 48.0f;

	static float r;
	static float k;
	static float l;
	static float wZero;

	public static float getRadius() {
		return r;
	}

	public static void setRadius(final float radius) {
		r = radius;
		final float h3 = r * r * r;
		k = K_FACTOR / (MathUtils.PI * h3);
		l = L_FACTOR / (MathUtils.PI * h3);
		wZero = W(new Vector3());
	}

	public static float W(final Vector3 rv) {
		final float rl = rv.len();
		final float q = rl / r;
		if (q <= 0.5f) {
			final float q2 = q * q;
			final float q3 = q2 * q;
			return k * (6.0f * q3 - 6.0f * q2 + 1.0f);
		} else {
			final float fac = 1.0f - q;
			return k * (2.0f * (fac * fac * fac));
		}
	}

	// FIXME (kotlin vector3 extension)
	public static Vector3 gradW(final Vector3 rv) {
		final float rl = rv.len();
		final float q = rl / r;
		if (rl > 1.0e-6) {
			final Vector3 gradq = rv.scl(1.0f / (rl * r));
			if (q <= 0.5f) {
				return gradq.scl(l * q * 3.0f * q - 2.0f);
			} else {
				final float fac = 1.0f - q;
				return gradq.scl(l * (-fac * fac));
			}
		}
		return Vector3.Zero;
	}

	public static float W_ZERO() { return wZero; }
}
