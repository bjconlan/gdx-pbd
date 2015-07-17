package pbd;

import com.badlogic.gdx.math.MathUtils;
import com.badlogic.gdx.math.Vector3;

/**
 * A non threadsafe sph kernel implementation (cubic kernel in cpp source)
 * as used by the position based dynamics class.
 *
 * @version ${project.version}
 * @since 1.1.0
 */
public class SPHKernel { // this should be a singleton
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
		float h3 = r * r * r;
		k = 8f / (MathUtils.PI * h3);
		l = 48f / (MathUtils.PI * h3);
		wZero = w(new Vector3());
	}

	public float wZero() {
		return wZero;
	}

	public float w(final Vector3 rv) {
		float rl = rv.len();
		float q = rl / r;
		if (q <= 0.5f) {
			float q2 = q * q;
			float q3 = q2 * q;
			return k * (6.0f * q3 - 6.0f * q2 + 1.0f);
		} else {
			float fac = 1.0f - q;
			return k * (2.0f * (fac * fac * fac));
		}
	}

	public Vector3 gradW(final Vector3 rv) {
		float rl = rv.len();
		float q = rl / r;
		if (rl > 1.0e-6) {
			Vector3 gradq = new Vector3(rv).scl(1.0f / (rl * r));
			if (q <= 0.5f) {
				return gradq.scl(l * q * (3.0f * q - 2.0f));
			} else {
				final float fac = 1.0f - q;
				return gradq.scl(l * (-fac * fac));
			}
		}
		return Vector3.Zero;
	}
}
