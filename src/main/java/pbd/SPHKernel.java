package pbd;

import com.badlogic.gdx.math.MathUtils;
import com.badlogic.gdx.math.Vector3;

public class SPHKernel {
	private float r;
	private float k;
	private float l;
	private float wZero;

	public SPHKernel(float radius) {
		setRadius(radius);
	}

	public float getRadius() {
		return r;
	}

	public void setRadius(final float radius) {
		r = radius;
		final float h3 = r * r * r;
		k = 8f / (MathUtils.PI * h3);
		l = 48f / (MathUtils.PI * h3);
		wZero = w(new Vector3());
	}

	public float w(final Vector3 rv) {
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
	public Vector3 gradW(final Vector3 rv) {
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

	public float wZero() { return wZero; }
}
