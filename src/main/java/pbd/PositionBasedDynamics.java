package pbd;

import com.badlogic.gdx.math.*;
import java.util.Arrays;

/**
 * @version ${project.version}
 * @since ${project.version}
 */
public class PositionBasedDynamics {
	private SPHKernel sphKernel;

	public PositionBasedDynamics(SPHKernel sphKernel) {
		this.sphKernel = sphKernel;
	}

	// helper read col of matrix
	static Vector3 col(final Matrix3 mtx, int col) {
		float[] m = mtx.val;
		switch (col) {
			case 0  : return new Vector3(m[Matrix3.M00], m[Matrix3.M10], m[Matrix3.M20]);
			case 1  : return new Vector3(m[Matrix3.M01], m[Matrix3.M11], m[Matrix3.M21]);
			default : return new Vector3(m[Matrix3.M02], m[Matrix3.M12], m[Matrix3.M22]);
		}
	}

	// helper set col
	static void setCol(final Vector3 v, Matrix3 mtx, int col) {
		float[] m = mtx.val;
		switch (col) {
			case 0  : m[Matrix3.M00] = v.x; m[Matrix3.M10] = v.y; m[Matrix3.M20] = v.z; return;
			case 1  : m[Matrix3.M01] = v.x; m[Matrix3.M11] = v.y; m[Matrix3.M21] = v.z; return;
			default : m[Matrix3.M02] = v.x; m[Matrix3.M12] = v.y; m[Matrix3.M22] = v.z; return;
		}
	}

	// helper read row of matrix
	static Vector3 row(final Matrix3 mtx, int row) {
		float[] m = mtx.val;
		switch (row) {
			case 0  : return new Vector3(m[Matrix3.M00], m[Matrix3.M01], m[Matrix3.M02]);
			case 1  : return new Vector3(m[Matrix3.M10], m[Matrix3.M11], m[Matrix3.M12]);
			default : return new Vector3(m[Matrix3.M20], m[Matrix3.M21], m[Matrix3.M22]);
		}
	}

	// helper set row
	static void setRow(final Vector3 v, Matrix3 mtx, int row) {
		float[] m = mtx.val;
		switch (row) {
			case 0  : m[Matrix3.M00] = v.x; m[Matrix3.M01] = v.y; m[Matrix3.M02] = v.z; return;
			case 1  : m[Matrix3.M10] = v.x; m[Matrix3.M11] = v.y; m[Matrix3.M12] = v.z; return;
			default : m[Matrix3.M20] = v.x; m[Matrix3.M21] = v.y; m[Matrix3.M22] = v.z; return;
		}
	}

	// NOTE we pretend this is a matrix2
	static float get(Affine2 m, int c, int r) {
		return c == 0
				? r == 0 ? m.m00 : m.m01
				: r == 0 ? m.m10 : m.m11;
	}

	static float get(Vector3 v, int c) {
		return c == 0 ? v.x : c == 1 ? v.y : v.z;
	}
	static void set(Vector3 v, int c, float val) {
		if (c == 0) { v.x = val; } else if (c == 1) { v.y = val; } else { v.z = val; }
	}

	// helper to get position in matrix
	static int mtx3Pos(int r, int c) {
		return (c * 3) + r;
	}

	static int mtx4Pos(int r, int c) {
		return (c * 4) * r;
	}

	static void jacobiRotate(Matrix3 a, Matrix3 r, int p, int q) {
		if (a.val[mtx3Pos(p, q)] == 0f) { // check this
			return;
		}

		float d = (a.val[mtx3Pos(p, p)] - a.val[mtx3Pos(q, q)]) / (2f * a.val[mtx3Pos(p, q)]);
		float t = 1.0f / (Math.abs(d) + (float)Math.sqrt(d * d + 1f));
		if (d < 0f) {
			t = -t;
		}
		float c = 1.0f / (float)Math.sqrt((t * t) + 1); // this can be moved up above the if
		float s = t * c;
		a.val[mtx3Pos(p, p)] += t * a.val[mtx3Pos(p, q)];
		a.val[mtx3Pos(q, q)] -= t * a.val[mtx3Pos(p, q)];
		a.val[mtx3Pos(p, q)] = a.val[mtx3Pos(q, p)] = 0f;

		for (int k = 0; k < 3; k++) {
			if (k != p && k != q) {
				float Akp = c * a.val[mtx3Pos(k, p)] + s * a.val[mtx3Pos(k, q)];
				float Akq = -s * a.val[mtx3Pos(k, p)] + c * a.val[mtx3Pos(k, q)];
				a.val[mtx3Pos(k, p)] = a.val[mtx3Pos(p, k)] = Akp;
				a.val[mtx3Pos(k, q)] = a.val[mtx3Pos(q, k)] = Akq;
			}
		}

		for (int k = 0; k < 3; k++) {
			float Rkp = c * r.val[mtx3Pos(k, p)] + s * r.val[mtx3Pos(k, q)];
			float Rkq = -s * r.val[mtx3Pos(k, p)] + c * r.val[mtx3Pos(k, q)];
			r.val[mtx3Pos(k, p)] = Rkp;
			r.val[mtx3Pos(k, q)] = Rkq;
		}
		// review with bullet btMatrix3x3 impl.
	}

	static void eigenDecomposition(final Matrix3 A, Matrix3 eigenVecs, Vector3 eigenVals) {
		final int numJacobiIterations = 10;
		final float epsilon = 1e-15f;

		Matrix3 D = new Matrix3(A);
		eigenVecs.idt();
		int iter = 0;
		while (iter < numJacobiIterations) { // iter < 10
			int p = 0;
			int q = 1;
			float max = Math.abs(D.val[Matrix3.M01]);
			float a = Math.abs(D.val[Matrix3.M02]);
			if (a > max) {
				p = 0; q = 2; max = a;
			}
			a = Math.abs(D.val[Matrix3.M12]);
			if (a > max) {
				p = 1; q = 2; max = a;
			}
			if (max < epsilon) {
				break;
			}
			jacobiRotate(A, eigenVecs, p, q);
			iter++;
		}
		eigenVals.x = D.val[Matrix3.M00];
		eigenVals.y = D.val[Matrix3.M11];
		eigenVals.z = D.val[Matrix3.M22];
	}

	static void polarDecomposition(final Matrix3 a, Matrix3 r, Matrix3 u, Matrix3 d) {
		Matrix3 AAT = new Matrix3();
		AAT.val[Matrix3.M00] = a.val[Matrix3.M00] * a.val[Matrix3.M00] + a.val[Matrix3.M01] * a.val[Matrix3.M01] + a.val[Matrix3.M02] * a.val[Matrix3.M02];
		AAT.val[Matrix3.M11] = a.val[Matrix3.M10] * a.val[Matrix3.M10] + a.val[Matrix3.M11] * a.val[Matrix3.M11] + a.val[Matrix3.M12] * a.val[Matrix3.M12];
		AAT.val[Matrix3.M22] = a.val[Matrix3.M20] * a.val[Matrix3.M20] + a.val[Matrix3.M21] * a.val[Matrix3.M21] + a.val[Matrix3.M22] * a.val[Matrix3.M22];

		AAT.val[Matrix3.M01] = a.val[Matrix3.M00] * a.val[Matrix3.M10] + a.val[Matrix3.M01] * a.val[Matrix3.M11] + a.val[Matrix3.M02] * a.val[Matrix3.M12];
		AAT.val[Matrix3.M02] = a.val[Matrix3.M00] * a.val[Matrix3.M20] + a.val[Matrix3.M01] * a.val[Matrix3.M21] + a.val[Matrix3.M02] * a.val[Matrix3.M22];
		AAT.val[Matrix3.M12] = a.val[Matrix3.M10] * a.val[Matrix3.M20] + a.val[Matrix3.M11] * a.val[Matrix3.M21] + a.val[Matrix3.M12] * a.val[Matrix3.M22];

		AAT.val[Matrix3.M10] = AAT.val[Matrix3.M01];
		AAT.val[Matrix3.M20] = AAT.val[Matrix3.M02];
		AAT.val[Matrix3.M21] = AAT.val[Matrix3.M12];

		r.idt();
		Vector3 eigenVals = new Vector3();
		eigenDecomposition(AAT, u, eigenVals);

		float d0 = (float)Math.sqrt(eigenVals.x);
		float d1 = (float)Math.sqrt(eigenVals.y);
		float d2 = (float)Math.sqrt(eigenVals.z);

		Arrays.fill(d.val, 0f);
		d.val[Matrix3.M00] = d0;
		d.val[Matrix3.M11] = d1;
		d.val[Matrix3.M22] = d2;

		final float eps = 1e-15f;

		float l0 = eigenVals.x <= eps ? 0f : 1.0f / d0;
		float l1 = eigenVals.y <= eps ? 0f : 1.0f / d1;
		float l2 = eigenVals.z <= eps ? 0f : 1.0f / d2;

		Matrix3 S1 = new Matrix3();
		S1.val[Matrix3.M00] = l0 * u.val[Matrix3.M00] * u.val[Matrix3.M00] + l1 * u.val[Matrix3.M01] * u.val[Matrix3.M01] + l2 * u.val[Matrix3.M02] * u.val[Matrix3.M02];
		S1.val[Matrix3.M11] = l0 * u.val[Matrix3.M10] * u.val[Matrix3.M10] + l1 * u.val[Matrix3.M11] * u.val[Matrix3.M11] + l2 * u.val[Matrix3.M12] * u.val[Matrix3.M12];
		S1.val[Matrix3.M22] = l0 * u.val[Matrix3.M20] * u.val[Matrix3.M20] + l1 * u.val[Matrix3.M21] * u.val[Matrix3.M21] + l2 * u.val[Matrix3.M22] * u.val[Matrix3.M22];

		S1.val[Matrix3.M01] = l0 * u.val[Matrix3.M00] * u.val[Matrix3.M10] + l1 * u.val[Matrix3.M01] * u.val[Matrix3.M11] + l2 * u.val[Matrix3.M02] * u.val[Matrix3.M12];
		S1.val[Matrix3.M02] = l0 * u.val[Matrix3.M00] * u.val[Matrix3.M20] + l1 * u.val[Matrix3.M01] * u.val[Matrix3.M21] + l2 * u.val[Matrix3.M02] * u.val[Matrix3.M22];
		S1.val[Matrix3.M12] = l0 * u.val[Matrix3.M10] * u.val[Matrix3.M20] + l1 * u.val[Matrix3.M11] * u.val[Matrix3.M21] + l2 * u.val[Matrix3.M12] * u.val[Matrix3.M22];

		S1.val[Matrix3.M10] = S1.val[Matrix3.M01];
		S1.val[Matrix3.M20] = S1.val[Matrix3.M02];
		S1.val[Matrix3.M21] = S1.val[Matrix3.M12];

		r = S1.mul(a);

		// col double check this... i thin kits wrong (rows, not colums)
		Vector3 c0 = col(r, 0);
		Vector3 c1 = col(r, 1);
		Vector3 c2 = col(r, 2);

		if (c0.len2() < eps) {
			c0 = c1.crs(c2);
		} else if (c1.len2() < eps) {
			c1 = new Vector3(c2).crs(c0);
		} else {
			c2 = new Vector3(c0).crs(c1);
		}

		setCol(c0, r, 0);
		setCol(c1, r, 1);
		setCol(c2, r, 2);
	}

	// find the largest sum for a col in the matrix (map > over sum of each col in mtx)
	static float oneNorm(final Matrix3 a) {
		final float sum1 = Math.abs(a.val[Matrix3.M00]) + Math.abs(a.val[Matrix3.M10]) + Math.abs(a.val[Matrix3.M20]);
		final float sum2 = Math.abs(a.val[Matrix3.M01]) + Math.abs(a.val[Matrix3.M11]) + Math.abs(a.val[Matrix3.M21]);
		final float sum3 = Math.abs(a.val[Matrix3.M02]) + Math.abs(a.val[Matrix3.M12]) + Math.abs(a.val[Matrix3.M22]);
		float maxSum = sum1;
		if (sum2 > maxSum) {
			maxSum = sum2;
		}
		if (sum3 > maxSum) {
			maxSum = sum3;
		}
		return maxSum;
	}

	// find the largest sum for a row in the matrix (map > over sum of each row in mtx)
	static float infNorm(final Matrix3 a) {
		final float sum1 = Math.abs(a.val[Matrix3.M00]) + Math.abs(a.val[Matrix3.M01]) + Math.abs(a.val[Matrix3.M02]);
		final float sum2 = Math.abs(a.val[Matrix3.M10]) + Math.abs(a.val[Matrix3.M11]) + Math.abs(a.val[Matrix3.M12]);
		final float sum3 = Math.abs(a.val[Matrix3.M20]) + Math.abs(a.val[Matrix3.M21]) + Math.abs(a.val[Matrix3.M22]);
		float maxSum = sum1;
		if (sum2 > maxSum) {
			maxSum = sum2;
		}
		if (sum3 > maxSum) {
			maxSum = sum3;
		}
		return maxSum;
	}

	public static void polarDecomposition2(final Matrix3 m, final float tolerance, Matrix3 r) {
		Matrix3 mT = m.transpose();
		float mOne = oneNorm(m);
		float mInf = infNorm(m);
		float eOne;
		Matrix3 mAdjTt = new Matrix3();
		Matrix3 eT = new Matrix3();

		do {
			setRow(row(mT, 1).crs(row(mT, 2)), mAdjTt, 0);
			setRow(row(mT, 2).crs(row(mT, 0)), mAdjTt, 1);
			setRow(row(mT, 0).crs(row(mT, 1)), mAdjTt, 2);

			float det = mT.val[Matrix3.M00] * mAdjTt.val[Matrix3.M00]
					+ mT.val[Matrix3.M01] * mAdjTt.val[Matrix3.M01]
					+ mT.val[Matrix3.M02] * mAdjTt.val[Matrix3.M02];

			if (Math.abs(det) < 1.0e-12) {
				int index = Integer.MAX_VALUE;
				for (int i = 0; i < 3; i++) {
					if (row(mAdjTt, i).len2() > 1.0e-12) {
						index = i;
						break;
					}
				}
				if (index == Integer.MAX_VALUE) {
					r.idt();
					return;
				} else {
					setRow(row(mT, (index + 1) % 3).crs(row(mT, (index + 2) % 3)), mT, index);
					setRow(row(mT, (index + 2) % 3).crs(row(mT, index)), mAdjTt, (index + 1) % 3);
					setRow(row(mT, index).crs(row(mT, (index + 1) % 3)), mAdjTt, (index + 2) % 3);
					Matrix3 m2 = new Matrix3(mT.val).transpose();
					mOne = oneNorm(m2);
					mInf = infNorm(m2);
					det = mT.val[Matrix3.M00] * mAdjTt.val[Matrix3.M00]
							+ mT.val[Matrix3.M01] * mAdjTt.val[Matrix3.M01]
							+ mT.val[Matrix3.M02] * mAdjTt.val[Matrix3.M02];
				}
			}
			float mAdjTOne = oneNorm(mAdjTt);
			float mAdjTInf = infNorm(mAdjTt);

			float gamma = (float)Math.sqrt((mAdjTOne * mAdjTInf) / (mOne * mInf)) / Math.abs(det);

			float g1 = gamma * 0.5f;
			float g2 = 0.5f / (gamma * det);

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					int ij = mtx3Pos(i, j);
					eT.val[ij] = mT.val[ij];
					mT.val[ij] = g1 * mT.val[ij] + g2 * mAdjTt.val[ij];
					eT.val[ij] -= mT.val[ij];
				}
			}

			eOne = oneNorm(eT);
			mOne = oneNorm(mT);
			mInf = oneNorm(mT);
		} while (eOne > mOne * tolerance);

		r.set(mT.transpose());
	}

	boolean solveDistanceConstraint(final Vector3 p0, final float invMass0,
	                                final Vector3 p1, final float invMass1,
	                                final float restLength,
	                                final float compressionStiffness,
	                                final float streatchStiffness,
	                                Vector3 corr0, Vector3 corr1) {
		float wSum = invMass0 + invMass1;
		if (wSum == 0f) {
			return false;
		}
		Vector3 n = new Vector3(p1).sub(p0);
		float d = n.len();
		n.nor();

		Vector3 corr = d < restLength
				? new Vector3(n).scl(compressionStiffness * (d - restLength) / wSum)
				: new Vector3(n).scl(streatchStiffness * (d - restLength) / wSum);
		return true;
	}

	boolean solveDihedralConstraint(final Vector3 p0, float invMass0,
	                                final Vector3 p1, float invMass1,
	                                final Vector3 p2, float invMass2,
	                                final Vector3 p3, float invMass3,
	                                float restAngle,
	                                float stiffness,
	                                Vector3 corr0, Vector3 corr1, Vector3 corr2, Vector3 corr3) {
		if (invMass0 == 0f && invMass1 == 0f) {
			return false;
		}

		Vector3 e = new Vector3(p3).sub(p2);
		float eLen = e.len();
		if (eLen < 1e-6f) {
			return false;
		}

		float invELen = 1.0f / eLen;
		Vector3 n1 = new Vector3(p2).sub(p0).crs(new Vector3(p3).sub(p0));
		n1.scl(1 / n1.len2());
		Vector3 n2 = new Vector3(p3).sub(p1).crs(new Vector3(p2).sub(p1));
		n2.scl(1 / n2.len2());

		Vector3 d0 = new Vector3(n1).scl(eLen);
		Vector3 d1 = new Vector3(n2).scl(eLen);
		Vector3 d2 = new Vector3(n1).scl(new Vector3(p0).sub(p3).dot(e) * invELen).add(
				new Vector3(n2).scl(new Vector3(p1).sub(p3).dot(e) * invELen));
		Vector3 d3 = new Vector3(n1).scl(new Vector3(p2).sub(p0).dot(e) * invELen).add(
				new Vector3(n2).scl(new Vector3(p2).sub(p1).dot(e) * invELen));

		n1.nor();
		n2.nor();
		float dot = n1.dot(n2); // inline

		if (dot < -1f) { dot = -1.0f; } // clamp
		if (dot > 1f) { dot = 1.0f; }
		float phi = (float)Math.acos(dot);

		float lambda =
				invMass0 * d0.len2() +
						invMass1 * d1.len2() +
				invMass2 * d2.len2() +
				invMass3 * d3.len2();

		if (lambda == 0f) {
			return false;
		}

		if (new Vector3(n1).crs(n2).dot(e) > 0f) {
			lambda = -lambda;
		}

		corr0 = d0.scl(- invMass0 * lambda); // REVIEW
		corr1 = d1.scl(- invMass1 * lambda);
		corr2 = d2.scl(- invMass2 * lambda);
		corr3 = d3.scl(- invMass3 * lambda);

		return true;
	}

	// FIXME remove lots of silly assignments which are never used
	boolean solveVolumeConstraint(final Vector3 p0, final float invMass0,
	                              final Vector3 p1, final float invMass1,
	                              final Vector3 p2, final float invMass2,
	                              final Vector3 p3, final float invMass3,
	                              final float restVolume,
	                              final float negVolumeStiffness,
	                              final float posVolumesStiffness,
	                              Vector3 corr0, Vector3 corr1, Vector3 corr2, Vector3 corr3) {
		Vector3 d1 = new Vector3(p1).sub(p0); // these can be removed and inlined FIXME
		Vector3 d2 = new Vector3(p2).sub(p0);
		Vector3 d3 = new Vector3(p3).sub(p0);
		float volume = 1.0f / 6.0f * (new Vector3(p1).sub(p0).crs(d2).dot(d3));

		corr0.setZero();
		corr1.setZero();
		corr2.setZero();
		corr3.setZero();

		if (posVolumesStiffness == 0.0f && volume > 0.0f) {
			return false;
		}
		if (negVolumeStiffness == 0.0f && volume < 0.0f) {
			return false;
		}

		Vector3 grad0 = new Vector3(p1).sub(p2).crs(new Vector3(p3).sub(p2));
		Vector3 grad1 = new Vector3(p2).sub(p0).crs(new Vector3(p3).sub(p0));
		Vector3 grad2 = new Vector3(p0).sub(p1).crs(new Vector3(p3).sub(p1));
		Vector3 grad3 = new Vector3(p1).sub(p0).crs(new Vector3(p2).sub(p0));

		float lambda = invMass0 * grad0.len2() +
				invMass1 * grad1.len2() +
				invMass2 * grad2.len2() +
				invMass3 * grad3.len2();

		if (Math.abs(lambda) < 1.0e-9) {
			return false;
		}

		lambda = volume < 0.0f
				? negVolumeStiffness * (volume - restVolume) / lambda
				: posVolumesStiffness * (volume - restVolume) / lambda;

		corr0 = grad0.scl(-lambda * invMass0); // we don't need to zero corr (we could just never use it
		corr1 = grad1.scl(-lambda * invMass1);
		corr2 = grad2.scl(-lambda * invMass2);
		corr3 = grad3.scl(-lambda * invMass3);

		return true;
	}

	float cotTheta(final Vector3 v, final Vector3 w) {
		final float cosTheta = v.dot(w);
		final float sinTheta = v.crs(w).len();
		return (cosTheta / sinTheta);
	}

	boolean computeQuadraticBendingMat(final Vector3 p0, final Vector3 p1, final Vector3 p2, final Vector3 p3,
	                                Matrix4 q) {
		Vector3 e0 = new Vector3(p3).sub(p0);
		Vector3 e1 = new Vector3(p2).sub(p0);
		Vector3 e2 = new Vector3(p3).sub(p0);
		Vector3 e3 = new Vector3(p2).sub(p1);
		Vector3 e4 = new Vector3(p3).sub(p1);

		Vector3 e0Neg = new Vector3(p0).scl(-1.0f);

		float c1 = cotTheta(e0, e1);
		float c2 = cotTheta(e0, e2);
		float c3 = cotTheta(e0Neg, e3);
		float c4 = cotTheta(e0Neg, e4);

		final float a0 = 0.5f * new Vector3(e0).crs(e1).len();
		final float a1 = 0.5f * new Vector3(e0).crs(e2).len();

		final float coef = -3f / (2f * (a0 + a1));
		final float[] k = { c3 + c4, c1 + c2, -c1 - c3, -c2 - c4 };
		final float[] k2 = { coef * k[0], coef * k[1], coef * k[2], coef * k[3] };

		float[] qv = q.val;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < i; j++) {
				qv[mtx4Pos(i, j)] = qv[mtx4Pos(j, i)] = k[i] * k2[j];
			}
			qv[mtx4Pos(i, i)] = k[i] * k2[i];
		}

		return true;
	}

	boolean solveIsometricBendingConstraint(final Vector3 p0, final float invMass0,
	                                        final Vector3 p1, final float invMass1,
	                                        final Vector3 p2, final float invMass2,
	                                        final Vector3 p3, final float invMass3,
	                                        final Matrix4 q,
	                                        final float stiffness,
	                                        Vector3 corr0, Vector3 corr1, Vector3 corr2, Vector3 corr3) {
		final Vector3[] x = { p2, p3, p0, p1 };
		float[] invMass = { invMass2, invMass3, invMass0, invMass1 };
		float energy = 0f;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				energy += q.val[mtx4Pos(i, j)] * x[j].dot(x[i]);
			}
		}
		energy *= 0.5f;

		Vector3[] gradC = { Vector3.Zero, Vector3.Zero, Vector3.Zero, Vector3.Zero };
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				gradC[j].add(new Vector3(x[i]).scl(q.val[mtx4Pos(j, i)]));
			}
		}

		float sumNormGradC = 0f;
		for (int i = 0; i < 4; i++) {
			if (invMass[i] != 0.0f) {
				sumNormGradC += invMass[i] * gradC[i].len2();
			}
		}

		if (Math.abs(sumNormGradC) > 1.0e-9) {
			final float s = energy / sumNormGradC;
			corr0.set(gradC[2].scl(-stiffness * s * invMass[2]));
			corr0.set(gradC[3].scl(-stiffness * s * invMass[3]));
			corr0.set(gradC[0].scl(-stiffness * s * invMass[0]));
			corr0.set(gradC[1].scl(-stiffness * s * invMass[1]));

			return true;
		}

		return false;
	}

	boolean solvedEdgePointDistConstraint(final Vector3 p, final float invMass,
	                                      final Vector3 p0, final float invMass0,
	                                      final Vector3 p1, final float invMass1,
	                                      final float restDist,
	                                      final float compressionStiffness,
	                                      final float stretchStiffness,
	                                      Vector3 corr, Vector3 corr0, Vector3 corr1) {
		float EPS = 1e-6f;
		Vector3 d = new Vector3(p1).sub(p0);
		float t;

		if (new Vector3(p0).sub(p1).len2() < (EPS * EPS)) {
			t = 0.5f;
		} else {
			float d2 = d.dot(d);
			t = d.dot(new Vector3(p).sub(p1).scl(1 / d2)); // optimize (inline an ?:)
			if (t < 0.0f) { // MathUtils.clamp(t, 0, 1.0f);
				t = 0.0f;
			} else if (t > 1.0f) {
				t = 1.0f;
			}
		}
		Vector3 q = new Vector3(d).scl(t).add(p0);
		Vector3 n = new Vector3(p).sub(q);
		float dist = n.len();
		n.nor();
		float c = dist - restDist;
		float b0 = 1.0f - t;
		float b1 = t;
		Vector3 grad = n;
		Vector3 grad0 = new Vector3(n).scl(-b0);
		Vector3 grad1 = new Vector3(n).scl(-b1);

		float s = invMass + invMass0 * b0 * b0 + invMass1 * b1 * b1;
		if (s == 0.0f) {
			return false;
		}
		s = c / s;
		s *= (c < 0.0f) ? compressionStiffness : stretchStiffness;

		if (s == 0.0f) {
			return false;
		}

		corr = grad.scl(-s * invMass);
		corr0 = grad0.scl(-s * invMass0);
		corr1 = grad1.scl(-s * invMass1);
		return true;
	}

	boolean solveTrianglePointDistConstraint(final Vector3 p, final float invMass,
	                                         final Vector3 p0, final float invMass0,
	                                         final Vector3 p1, final float invMass1,
	                                         final Vector3 p2, final float invMass2,
	                                         final float restDist,
	                                         final float compressionStiffness,
	                                         final float stretchStiffness,
	                                         Vector3 corr, Vector3 corr0, Vector3 corr1, Vector3 corr2) {
		float b0 = 1.0f / 3.0f;
		float b1 = b0;
		float b2 = b0;

		Vector3 d1 = new Vector3(p1).sub(p0);
		Vector3 d2 = new Vector3(p2).sub(p0);
		Vector3 pp0 = new Vector3(p).sub(p0);
		float a = d1.dot(d1);
		float b = d2.dot(d1);
		float c = pp0.dot(d1);
		float d = b;
		float e = d2.dot(d2);
		float f = pp0.dot(d2);
		float det = a * e - b * d;

		if (det != 0.0f) {
			float s = (c * e - b * f) / det;
			float t = (a * f - c * d) / det;
			b0 = 1.0f - s - t;
			b1 = s;
			b2 = t;
			if (b0 < 0.0f) { // on edge 1-2 (optimize lot of common init code)
				Vector3 dv = new Vector3(p2).sub(p1);
				float d2s = dv.dot(dv);
				float ts = (d2s == 0.0f) ? 0.5f : dv.dot(new Vector3(p).sub(p1)) / d2s;
				if (ts < 0.0f) {
					ts = 0.0f;
				} // clamp
				if (ts > 1.0f) {
					ts = 1.0f;
				}
				b0 = 0f;
				b1 = 1.0f - ts;
				b2 = ts;
			} else if (b1 < 0.0f) { // on edge 2-0
				Vector3 dv = new Vector3(p0).sub(p2);
				float d2s = dv.dot(dv);
				float ts = (d2s == 0.0f) ? 0.5f : dv.dot(new Vector3(p).sub(p2)) / d2s;
				if (ts < 0.0f) {
					ts = 0.0f;
				} // clamp
				if (ts > 1.0f) {
					ts = 1.0f;
				}
				b1 = 0.0f;
				b2 = (1.0f - ts);
				b0 = ts;
			} else if (b2 < 0.0f) { // on edge 0-1
				Vector3 dv = new Vector3(p1).sub(p0);
				float d2s = dv.dot(dv);
				float ts = (d2s == 0.0f) ? 0.5f : dv.dot(new Vector3(p).sub(p0)) / d2s;
				if (ts < 0.0f) {
					ts = 0.0f;
				} // clamp
				if (ts > 1.0f) {
					ts = 1.0f;
				}
				b1 = 0.0f;
				b2 = (1.0f - ts);
				b0 = ts;
			}
		}
		Vector3 q = new Vector3(p0).scl(b0).add( // optimize?
				new Vector3(p1).scl(b1)).add(
				new Vector3(p2).scl(p2));
		Vector3 n = new Vector3(p).sub(q);
		float dist = n.len();
		n.nor();
		float C = dist - restDist; //Validate this
		Vector3 grad = new Vector3(n);
		Vector3 grad0 = new Vector3(n).scl(-b0);
		Vector3 grad1 = new Vector3(n).scl(-b1);
		Vector3 grad2 = new Vector3(n).scl(-b2);

		float s = invMass
				+ invMass0 * b0 * b0
				+ invMass1 * b1 * b1
				+ invMass2 * b2 * b2;
		if (s == 0f) {
			return false;
		}
		s = C / s;
		if (C < 0.0f) {
			s *= compressionStiffness;
		} else {
			s *= stretchStiffness;
		}
		if (s == 0.0f) {
			return false;
		}

		corr = grad.scl(-s * invMass);
		corr0 = grad.scl(-s * invMass0);
		corr1 = grad.scl(-s * invMass1);
		corr2 = grad.scl(-s * invMass2);
		return true;
	}

	boolean solveEdgeEdgeDistConstraint(final Vector3 p0, final float invMass0,
	                                    final Vector3 p1, final float invMass1,
	                                    final Vector3 p2, final float invMass2,
	                                    final Vector3 p3, final float invMass3,
	                                    final float restDist,
	                                    final float compressionStiffness,
	                                    final float stretchStiffness,
	                                    Vector3 corr0, Vector3 corr1, Vector3 corr2, Vector3 corr3) {
		Vector3 d0 = new Vector3(p1).sub(p0);
		Vector3 d1 = new Vector3(p3).sub(p2);

		float a = d0.len2();
		float b = -d0.dot(d0);
		float c = d0.dot(d0);
		float d = -d1.len2();
		float e = new Vector3(p2).sub(p0).dot(d0);
		float f = new Vector3(p2).sub(p0).dot(d1);
		float det = a * d - b * c;
		float s, t;
		if (det != 0.0f) {
			det = 1.0f / det;
			s = (e * d - b * f) * det;
			t = (a * f - e * c) * det;
		} else {
			float s0 = p0.dot(d0);
			float s1 = p1.dot(d0);
			float t0 = p2.dot(p0);
			float t1 = p3.dot(p0);
			boolean flip0 = false;
			boolean flip1 = false;

			if (s0 > s1) { float swp = s0; s0 = s1; s1 = swp; flip0 = true; }
			if (t0 > t1) { float swp = t0; t0 = t1; t1 = swp; flip1 = true; }

			if (s0 >= t1) {
				s = !flip0 ? 0.0f : 1.0f;
				t = !flip1 ? 1.0f : 0.0f;
			} else if (t0 >= s1) {
				s = !flip0 ? 1.0f : 0.0f;
				t = !flip1 ? 0.0f : 1.0f;
			} else {
				float mid = (s0 > t0) ? (s0 + t1) * 0.5f : (t0 + s1) * 0.5f;
				s = (s0 == s1) ? 0.5f : (mid - s0) / (s1 - s0);
				t = (t0 == t1) ? 0.5f : (mid - t0) / (t1 - t0);
			}
		}

		if (s < 0.0f) { s = 0.0f; }
		if (s > 1.0f) { s = 1.0f; }
		if (t < 0.0f) { t = 0.0f; }
		if (t > 1.0f) { t = 1.0f; }

		float b0 = 1.0f - s;
		float b1 = s;
		float b2 = 1.0f - t;
		float b3 = t;

		Vector3 q0 = new Vector3(p0).scl(b0).add(new Vector3(p1).scl(b1));
		Vector3 q1 = new Vector3(p2).scl(b2).add(new Vector3(p3).scl(b3));
		Vector3 n = new Vector3(q0).sub(q1);
		float dist = n.len();
		n.nor();
		float C = dist - restDist;
		Vector3 grad0 = new Vector3(n).scl(b0);
		Vector3 grad1 = new Vector3(n).scl(b1);
		Vector3 grad2 = new Vector3(n).scl(-b2);
		Vector3 grad3 = new Vector3(n).scl(-b3);

		s = invMass0 * b0 * b0 + invMass1 * b1 * b1 + invMass2 * b2 * b2 + invMass3 * b3 * b3;
		if (s == 0.0f) {
			return false;
		}

		s = C / s;
		s *= (C < 0.0f) ? compressionStiffness : stretchStiffness;

		if (s == 0.0f) {
			return false;
		}

		corr0 = grad0.scl(-s * invMass0);
		corr1 = grad1.scl(-s * invMass1);
		corr2 = grad2.scl(-s * invMass2);
		corr3 = grad3.scl(-s * invMass3);
		return true;
	}

	boolean computeShapeMatchingRestInfo(final Vector3[] x0, final float invMasses[], int numPoints,
	                                     Vector3 restCm, Matrix3 invRestMat) {
		final float eps = 1e-6f;
		invRestMat.idt();

		// center of mass
		restCm.setZero();
		float wSum = 0.0f;
		for (int i = 0; i < numPoints; i++) {
			float wi = 1.0f / (invMasses[i] + eps);
			restCm.add(new Vector3(x0[i]).scl(wi));
			wSum += wi;
		}
		if (wSum == 0.0f) {
			return false;
		}
		restCm.scl(1 / wSum);

		// A
		Matrix3 A = new Matrix3(new float[9]);
		for (int i = 0; i < numPoints; i++) {
			final Vector3 qi = new Vector3(x0[i]).sub(restCm);
			float wi = 1.0f / (invMasses[i] + eps);
			float x2 = wi * qi.x * qi.x;
			float y2 = wi * qi.y * qi.y;
			float z2 = wi * qi.z * qi.z;
			float xy = wi * qi.x * qi.y;
			float xz = wi * qi.x * qi.z;
			float yz = wi * qi.y * qi.z;
			float[] a = A.val;
			a[Matrix3.M00] += x2; a[Matrix3.M01] += xy; a[Matrix3.M02] += xz;
			a[Matrix3.M10] += xy; a[Matrix3.M11] += y2; a[Matrix3.M12] += yz;
			a[Matrix3.M20] += xz; a[Matrix3.M21] += yz; a[Matrix3.M22] += z2;
		}
		float det = A.det();
		if (Math.abs(det) > 1.0e-9) {
			invRestMat = A.inv();
			return true;
		}

		return false;
	}

	boolean solveShapeMatchingConstraint(final Vector3[] x0, final Vector3[] x, final float[] invMasses, final int numPoints,
	                                     final Vector3 restCm,
	                                     final Matrix3 invRestMat,
	                                     final float stiffness,
	                                     final boolean allowStretch,
	                                     Vector3[] corr, Matrix3 rot) {
		final float eps = 1e-6f;
		for (Vector3 c : corr) {
			c.setZero();
		}

		Vector3 cm = new Vector3();
		float wSum = 0.0f;
		for (int i = 0; i < invMasses.length; i++) {
			float wi = 1.0f / (invMasses[i] + eps);
			cm.add(new Vector3(x[i]).scl(wi));
			wSum += wi;
		}
		if (wSum == 0.0f) {
			return false;
		}
		cm.scl(1 / wSum);

		// A
		Matrix3 mat = new Matrix3(new float[9]);
		for (int i = 0; i < invMasses.length; i++) {
			Vector3 q = new Vector3(x0[i]).sub(restCm);
			Vector3 p = new Vector3(x[i]).sub(cm);

			float w = 1.0f / (invMasses[i] + eps);
			p.scl(w);

			float[] m = mat.val;
			m[Matrix3.M00] += p.x * q.x; m[Matrix3.M01] += p.x * q.y; m[Matrix3.M02] += p.x * q.z;
			m[Matrix3.M10] += p.y * q.x; m[Matrix3.M11] += p.y * q.y; m[Matrix3.M12] += p.y * q.z;
			m[Matrix3.M20] += p.z * q.x; m[Matrix3.M21] += p.z * q.y; m[Matrix3.M22] += p.z * q.z;
		}

		mat.mul(invRestMat);

		// REVIEW
		Matrix3 R, U, D;
		R = mat;
		if (!allowStretch) {
			polarDecomposition2(mat, 1e-6f, R);
		}
		for (int i = 0; i < corr.length; i++) {
			Vector3 goal = new Vector3(x0[i]).sub(restCm).mul(R).add(cm);
			corr[i] = goal.sub(x[i]).scl(stiffness);
		}

		if (rot != null) {
			rot = R;
		}

		return true;
	}

	boolean computeStrainTriangleInvRestMat(final Vector3 p0, final Vector3 p1, final Vector3 p2,
	                                        Matrix2 invRestMat) {
		float a = p1.x - p0.x; float b = p2.x - p0.x;
		float c = p1.y - p0.y; float d = p2.y - p0.y;

		float det = a * d - b * c;
		if (Math.abs(det) < 1.0e-9) {
			return false;
		}

		float s = 1.0f / det;
		invRestMat.val[Matrix2.M00] = d * s; invRestMat.val[Matrix2.M01] = -b * s;
		invRestMat.val[Matrix2.M10] = -c * s; invRestMat.val[Matrix2.M11] = a * s;

		return true;
	}

	boolean solveStrainTriangleConstraint(final Vector3 p0, final float invMass0,
	                                      final Vector3 p1, final float invMass1,
	                                      final Vector3 p2, final float invMass2,
	                                      final Matrix3 invRestMat, // 2d Matrix
	                                      final float xxStiffness,
	                                      final float yyStiffness,
	                                      final float xyStiffness,
	                                      final boolean normalizeStretch,
	                                      final boolean normalizeShear,
	                                      Vector3 corr0, Vector3 corr1, Vector3 corr2) {
		Vector3[] c = new Vector3[3];
		c[0] = new Vector3(invRestMat.val[Matrix3.M00], invRestMat.val[Matrix3.M10], 0.0f);
		c[1] = new Vector3(invRestMat.val[Matrix3.M01], invRestMat.val[Matrix3.M11], 0.0f);

		Vector3[] r = new Vector3[3];

		corr0.setZero();
		corr1.setZero();
		corr2.setZero();

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j <= i; j++) {
				// alternative is Jacobi

				//Gauss - Seidel
				r[0] = new Vector3((p1.x + corr1.x) - (p0.x + corr0.x), (p2.x + corr2.x) - (p0.x + corr0.x), 0f);
				r[1] = new Vector3((p1.y + corr1.y) - (p0.y + corr0.y), (p2.y + corr2.y) - (p0.y + corr0.y), 0f);
				r[2] = new Vector3((p1.z + corr1.z) - (p0.z + corr0.z), (p2.z + corr2.z) - (p0.z + corr0.z), 0f);


				float Sij = 0.0f;
				for (int k = 0; k < 3; k++) {
					Sij += r[k].dot(c[i]) * r[k].dot(c[j]);
				}

				Vector3[] d = new Vector3[3];
				d[0] = Vector3.Zero;

				for (int k = 0; k < 2; k++) {
					d[k + 1] = new Vector3(r[0].dot(c[j]), r[1].dot(c[j]), r[2].dot(c[j])).scl(invRestMat.val[mtx3Pos(k, i)]);
					d[k + 1] = new Vector3(r[0].dot(c[i]), r[1].dot(c[i]), r[2].dot(c[i])).scl(invRestMat.val[mtx3Pos(k, j)]);
					d[0].sub(d[k + 1]);
				}

				if (i != j && normalizeShear) {
					float fi2 = 0.0f;
					float fj2 = 0.0f;
					for (int k = 0; k < 3; k++) {
						fi2 += r[k].dot(c[i]) * r[k].dot(c[i]);
						fj2 += r[k].dot(c[k]) * r[k].dot(c[j]);
					}
					float fi = (float) Math.sqrt(fi2);
					float fj = (float) Math.sqrt(fj2);

					d[0] = Vector3.Zero;
					float s = Sij / (fi2 * fi * fj2 * fj);
					for (int k = 0; k < 2; k++) {
						d[k + 1].scl(1 / (fi * fj));
						d[k + 1].sub(new Vector3(r[0].dot(c[i]), r[1].dot(c[i]), r[2].dot(c[i])).scl(invRestMat.val[mtx3Pos(k, i)] * s * fj * fj));
						d[k + 1].sub(new Vector3(r[0].dot(c[j]), r[1].dot(c[j]), r[2].dot(c[j])).scl(invRestMat.val[mtx3Pos(k, j)] * s * fi * fi));
						d[0].sub(d[k + 1]);
					}
					Sij = Sij / (fi * fj);
				}
				float lambda = invMass0 * d[0].len2()
						+ invMass1 * d[1].len2()
						+ invMass2 * d[2].len2();

				if (lambda == 0.0f) {
					continue;
				}

				if (i == 0 && j == 0) {
					if (normalizeStretch) {
						float s = (float) Math.sqrt(Sij);
						lambda = 2.0f * s * (s - 1.0f) / lambda * xxStiffness;
					} else {
						lambda = (Sij - 1.0f) / lambda * xxStiffness;
					}
				} else if (i == 1 && j == 1) {
					if (normalizeStretch) {
						float s = (float) Math.sqrt(Sij);
						lambda = 2.0f * s * (s - 1.0f) / lambda * yyStiffness;
					} else {
						lambda = (Sij - 1.0f) / lambda * yyStiffness;
					}
				} else {
					lambda = Sij / lambda * xyStiffness;
				}

				corr0.sub(d[0].scl(lambda * invMass0));
				corr1.sub(d[1].scl(lambda * invMass1));
				corr2.sub(d[2].scl(lambda * invMass2));
			}
		}
		return true;
	}

	boolean computeStrainTetraInvRestMat(final Vector3 p0, final Vector3 p1, final Vector3 p2, final Vector3 p3,
	                                     Matrix3 invRestMat) {
		Matrix3 m = new Matrix3();
		setCol(new Vector3(p1).sub(p0), m, 0);
		setCol(new Vector3(p2).sub(p1), m, 1);
		setCol(new Vector3(p3).sub(p0), m, 2);

		float det = m.det();
		if (Math.abs(det) > 1.0e-9) {
			invRestMat = m.inv();
			return true;
		}
		return false;
	}

	boolean solveStraintTetraConstraint(final Vector3 p0, final float invMass0,
	                                    final Vector3 p1, final float invMass1,
	                                    final Vector3 p2, final float invMass2,
	                                    final Vector3 p3, final float invMass3,
	                                    final Matrix3 invRestMat,
	                                    final Vector3 stretchStiffness,
	                                    final Vector3 shearStiffness,
	                                    final boolean normalizeStretch,
	                                    final boolean normalizeShear,
	                                    Vector3 corr0, Vector3 corr1, Vector3 corr2, Vector3 corr3) {
		corr0.setZero();
		corr1.setZero();
		corr2.setZero();
		corr3.setZero();

		Vector3[] c = new Vector3[] { col(invRestMat, 0), col(invRestMat, 1), col(invRestMat, 2) };

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j <= i; j++) {
				Matrix3 p = new Matrix3();
				// Gauss -Seidel (alternative is jacobi)
				Vector3 delta = new Vector3(p0).add(corr0);
				setCol(new Vector3(p1).add(corr1).sub(delta), p, 0);
				setCol(new Vector3(p2).add(corr2).sub(delta), p, 1);
				setCol(new Vector3(p3).add(corr3).sub(delta), p, 2);

				Vector3 fi = new Vector3(c[i]).mul(p);
				Vector3 fj = new Vector3(c[j]).mul(p);

				float Sij = fi.dot(fj);
				float wi = 0, wj = 0, s1= 0, s3 = 0;
				if (normalizeShear && i != j) {
					wi = fi.len(); // isn't this a well known algo? (1/vA.len * vV.len)^3
					wj = fj.len();
					s1 = 1.0f / (wi * wj);
					s3 = s1 * s1 * s1;
				}

				Vector3[] d = new Vector3[4];
				d[0] = Vector3.Zero;
				for (int k = 0; k < 3; k++) {
					d[k+1] = new Vector3(fj).scl(invRestMat.val[mtx3Pos(k, i)]).add(new Vector3(fi).scl(invRestMat.val[mtx3Pos(k, j)]));

					if (normalizeShear && i != j) {
						d[k+1] = new Vector3(d[k+1]).scl(s1).sub(Sij * s3)
								.scl(new Vector3(fi).scl(wj * wj * invRestMat.val[mtx3Pos(k, i)])
								.add(new Vector3(fj).scl(wi * wi * invRestMat.val[mtx3Pos(k, j)])));
					}
					d[0].sub(d[k+1]);
				}

				if (normalizeShear && i != j) {
					Sij *= s1;
				}

				float lambda = invMass0 * d[0].len2()
						+ invMass1 * d[1].len2()
						+ invMass2 * d[2].len2()
						+ invMass3 * d[3].len2();

				if (Math.abs(lambda) < 1e-6f) {
					continue;
				}

				if (i == j) {
					float ss = i == 0 ? stretchStiffness.x : i == 1 ? stretchStiffness.y : stretchStiffness.z;
					if (normalizeStretch) {
						float s = (float) Math.sqrt(Sij);
						lambda = 2.0f * s * (s - 1.0f) / lambda * ss;
					} else {
						lambda = (Sij - 1.0f) / lambda * ss;
					}
				} else {
					int idx = i + j - 1;
					float ss = idx == 0 ? shearStiffness.x : idx == 1 ? shearStiffness.y : shearStiffness.z;
					lambda = Sij / lambda * ss;
				}
				corr0.sub(d[0].scl(lambda * invMass0));
				corr1.sub(d[1].scl(lambda * invMass1));
				corr2.sub(d[2].scl(lambda * invMass2));
				corr2.sub(d[3].scl(lambda * invMass3));
			}
		}

		return true;
	}

	boolean computeFEMTriangleInvRestMat(final Vector3 p0, final Vector3 p1, final Vector3 p2,
	                                     float area, Matrix2 invRestMat) { // FIXME ref
		Vector3 normal0 = new Vector3(p1).sub(p0).crs(new Vector3(p2).sub(p0));
		area = normal0.len() * 0.0f; // FIXME no prim refs
		Vector3 axis0_1 = new Vector3(p1).sub(p0).nor();
		Vector3 axis0_2 = new Vector3(normal0).crs(axis0_1).nor();

		Vector2[] p = new Vector2[3];
		p[0] = new Vector2(p0.dot(axis0_2), p0.dot(axis0_1));
		p[1] = new Vector2(p1.dot(axis0_2), p1.dot(axis0_1));
		p[2] = new Vector2(p2.dot(axis0_2), p2.dot(axis0_1));

		Matrix2 pMtx = new Matrix2();
		pMtx.val[Matrix2.M00] = p[0].x - p[2].x;
		pMtx.val[Matrix2.M10] = p[0].y - p[2].y;
		pMtx.val[Matrix2.M01] = p[1].x - p[2].x;
		pMtx.val[Matrix2.M11] = p[1].y - p[2].y;

		final float det = pMtx.det();
		if (Math.abs(det) > 1.0e-9) {
			invRestMat = pMtx.inv(); // DOUBLE CHECK!
			return true;
		}

		return false;
	}

	boolean solveFEMTriangleConstraint(final Vector3 p0, final float invMass0,
	                                   final Vector3 p1, final float invMass1,
	                                   final Vector3 p2, final float invMass2,
	                                   final float area,
	                                   final Matrix2 invRestMat,
	                                   final float youngsModulusX,
	                                   final float youngsModulusY,
	                                   final float youngsModulusShear,
	                                   final float poissonRatioXY,
	                                   final float poissonRatioYX,
	                                   Vector3 corr0, Vector3 corr1, Vector3 corr2) {
		Matrix3 c = new Matrix3(new float[9]);
		c.val[Matrix3.M00] = youngsModulusX / (1.0f - poissonRatioXY * poissonRatioYX);
		c.val[Matrix3.M01] = youngsModulusX * poissonRatioYX / (1.0f - poissonRatioXY * poissonRatioYX);
		c.val[Matrix3.M11] = youngsModulusY / (1.0f - poissonRatioXY * poissonRatioYX);
		c.val[Matrix3.M10] = youngsModulusY * poissonRatioXY / (1.0f - poissonRatioXY * poissonRatioYX);
		c.val[Matrix3.M22] = youngsModulusShear;

		// determin partial x/partial m_i
		Affine2 f = new Affine2(); // mtx 3*2; REVIEW
		final Vector3 p13 = new Vector3(p0).sub(p2);
		final Vector3 p23 = new Vector3(p1).sub(p2);
		f.m00 = p13.x * invRestMat.val[Matrix2.M00] + p23.x * invRestMat.val[Matrix2.M10];
		f.m01 = p13.x * invRestMat.val[Matrix2.M01] + p23.x * invRestMat.val[Matrix2.M11];
		f.m10 = p13.y * invRestMat.val[Matrix2.M00] + p23.y * invRestMat.val[Matrix2.M10];
		f.m11 = p13.y * invRestMat.val[Matrix2.M01] + p23.y * invRestMat.val[Matrix2.M11];
		f.m02 = p13.z * invRestMat.val[Matrix2.M00] + p23.z * invRestMat.val[Matrix2.M10]; // review (col major)
		f.m12 = p13.z * invRestMat.val[Matrix2.M01] + p23.z * invRestMat.val[Matrix2.M11]; // review (col major)

		// epsilon = 0.5(F^T * F - I)
		Matrix2 epsilon = new Matrix2(); // Matrix2 REVIEW
		epsilon.val[Matrix2.M00] = 0.5f * (f.m00  * f.m00 + f.m10 * f.m10 + f.m02 * f.m02 - 1.0f);
		epsilon.val[Matrix2.M11] = 0.5f * (f.m01  * f.m01 + f.m11 * f.m11 + f.m12 * f.m12 - 1.0f);
		epsilon.val[Matrix2.M01] = epsilon.val[Matrix2.M10] = 0.5f * (f.m00  * f.m01 + f.m10 * f.m11 + f.m02 * f.m12);

		// p(f) = det(f) * c * e * f^-t => E = green strain
		Matrix2 stress = new Matrix2(); // Matrix2 Review
		stress.val[Matrix2.M00] = c.val[Matrix3.M00] * epsilon.val[Matrix2.M00] + c.val[Matrix3.M01] * epsilon.val[Matrix2.M11] + c.val[Matrix3.M02] * epsilon.val[Matrix2.M01];
		stress.val[Matrix2.M11] = c.val[Matrix3.M10] * epsilon.val[Matrix2.M00] + c.val[Matrix3.M11] * epsilon.val[Matrix2.M11] + c.val[Matrix3.M12] * epsilon.val[Matrix2.M01];
		stress.val[Matrix2.M01] = stress.val[Matrix2.M10] = c.val[Matrix3.M20] * epsilon.val[Matrix2.M00] + c.val[Matrix3.M21] * epsilon.val[Matrix2.M11] + c.val[Matrix3.M22] * epsilon.val[Matrix2.M01];

		final Affine2 piolaKirchhoffStres = new Affine2(f).mul(stress.toAffine2()); // Matrix2 Review
		float psi = 0.0f;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; i++) {
				psi += epsilon.get(i, j) * stress.get(i, j);
			}
		}
		psi = 0.5f * psi;
		float energy = area * psi;

		// compute gradient FIXME mul & transpose of affine are not symetrical
		Affine2 h = new Affine2(piolaKirchhoffStres).mul(invRestMat.transpose().toAffine2()); // review
		Vector3[] gradC = new Vector3[3];
		for (int i = 0; i < 3; i++) {
			float a0 = get(h, i, 0);
			float a1 = get(h, i, 0);
			switch (i) {
				case 0:
					gradC[0].x = a0;
					gradC[1].x = a1;
					break;
				case 1:
					gradC[0].y = a0;
					gradC[1].y = a1;
					break;
				default:
					gradC[0].z = a0;
					gradC[0].z = a1;
			}
		}
		gradC[2] = new Vector3(gradC[0]).scl(-1).sub(gradC[1]);

		float sumNormGradC = invMass0 * gradC[0].len2();
		sumNormGradC += invMass1 * gradC[1].len2();
		sumNormGradC += invMass2 * gradC[2].len2();

		if (Math.abs(sumNormGradC) > 1.0e-9) {
			float s = energy / sumNormGradC;
			corr0 = gradC[0].scl(-s * invMass0); // review (does - go on both?)
			corr1 = gradC[1].scl(-s * invMass1);
			corr2 = gradC[2].scl(-s * invMass2);
			return true;
		}
		return false;
	}

	boolean computeFEMTetraInvRestMat(final Vector3 p0, final Vector3 p1, final Vector3 p2, final Vector3 p3,
	                                  float volume, Matrix3 invRestMatrix) {
		Vector3 p3_p1 = new Vector3(p3).sub(p0);
		volume = Math.abs((1.0f / 6.0f) * new Vector3(p3).sub(p0).dot(new Vector3(p2).sub(p0).crs(new Vector3(p1).sub(p0))));

		Matrix3 m = new Matrix3();
		setCol(new Vector3(p0).sub(p3), m, 0);
		setCol(new Vector3(p1).sub(p3), m, 1);
		setCol(new Vector3(p2).sub(p3), m, 2);

		float det = m.det();
		if (Math.abs(det) > 1.0e-9) {
			invRestMatrix = m.inv();
			return true;
		}

		return false;
	}

	static void computeGreenStrainAndPiolaStress(
	        final Vector3 x1, final Vector3 x2, final Vector3 x3, final Vector3 x4,
	        final Matrix3 invRestMat, final float restVolume, final float mu, final float lambda,
	        Matrix3 epsilon, Matrix3 sigma, float energy) {
		Matrix3 f = new Matrix3();
		final Vector3 p14 = new Vector3(x1).sub(x4);
		final Vector3 p24 = new Vector3(x2).sub(x4);
		final Vector3 p34 = new Vector3(x3).sub(x4);

		float[] fv = f.val;
		float[] irm = invRestMat.val;
		fv[Matrix3.M00] = p14.x * irm[Matrix3.M00] + p24.x * irm[Matrix3.M10] + p34.x * irm[Matrix4.M20];
		fv[Matrix3.M01] = p14.x * irm[Matrix3.M01] + p24.x * irm[Matrix3.M11] + p34.x * irm[Matrix4.M21];
		fv[Matrix3.M02] = p14.x * irm[Matrix3.M02] + p24.x * irm[Matrix3.M12] + p34.x * irm[Matrix4.M22];
		fv[Matrix3.M10] = p14.y * irm[Matrix3.M00] + p24.y * irm[Matrix3.M10] + p34.y * irm[Matrix4.M20];
		fv[Matrix3.M11] = p14.y * irm[Matrix3.M01] + p24.y * irm[Matrix3.M11] + p34.y * irm[Matrix4.M21];
		fv[Matrix3.M12] = p14.y * irm[Matrix3.M02] + p24.y * irm[Matrix3.M12] + p34.y * irm[Matrix4.M22];
		fv[Matrix3.M20] = p14.z * irm[Matrix3.M00] + p24.z * irm[Matrix3.M10] + p34.z * irm[Matrix4.M20];
		fv[Matrix3.M21] = p14.z * irm[Matrix3.M01] + p24.z * irm[Matrix3.M11] + p34.z * irm[Matrix4.M21];
		fv[Matrix3.M22] = p14.z * irm[Matrix3.M02] + p24.z * irm[Matrix3.M12] + p34.z * irm[Matrix4.M22];

		// epsilon = 1/2 F^T F - I
		epsilon.val[Matrix3.M00] = 0.5f * (fv[Matrix3.M00] * fv[Matrix3.M00] + fv[Matrix3.M10] * fv[Matrix3.M10] + fv[Matrix3.M20] * fv[Matrix3.M20] - 1.0f); //xx
		epsilon.val[Matrix3.M11] = 0.5f * (fv[Matrix3.M01] * fv[Matrix3.M01] + fv[Matrix3.M11] * fv[Matrix3.M11] + fv[Matrix3.M21] * fv[Matrix3.M21] - 1.0f); //yy
		epsilon.val[Matrix3.M22] = 0.5f * (fv[Matrix3.M02] * fv[Matrix3.M02] + fv[Matrix3.M12] * fv[Matrix3.M12] + fv[Matrix3.M22] * fv[Matrix3.M22] - 1.0f); //zz
		epsilon.val[Matrix3.M01] = 0.5f * (fv[Matrix3.M00] * fv[Matrix3.M01] + fv[Matrix3.M10] * fv[Matrix3.M11] + fv[Matrix3.M20] * fv[Matrix3.M21]); // xy
		epsilon.val[Matrix3.M02] = 0.5f * (fv[Matrix3.M00] * fv[Matrix3.M02] + fv[Matrix3.M10] * fv[Matrix3.M12] + fv[Matrix3.M20] * fv[Matrix3.M22]); // xz
		epsilon.val[Matrix3.M12] = 0.5f * (fv[Matrix3.M01] * fv[Matrix3.M02] + fv[Matrix3.M11] * fv[Matrix3.M12] + fv[Matrix3.M21] * fv[Matrix3.M22]); // yz
		epsilon.val[Matrix3.M10] = epsilon.val[Matrix3.M01];
		epsilon.val[Matrix3.M20] = epsilon.val[Matrix3.M02];
		epsilon.val[Matrix3.M21] = epsilon.val[Matrix4.M12];

		final float trace = epsilon.val[Matrix3.M00] + epsilon.val[Matrix3.M11] + epsilon.val[Matrix3.M22];
		final float ltrace = lambda * trace;
		sigma = epsilon.scl(2.0f * mu);
		sigma.val[Matrix3.M00] += ltrace;
		sigma.val[Matrix3.M11] += ltrace;
		sigma.val[Matrix3.M22] += ltrace;
		sigma = f.mul(sigma);

		float psi = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				psi += epsilon.val[mtx3Pos(i, j)] * epsilon.val[mtx3Pos(i, j)];
			}
		}
		psi = mu * psi + 0.5f * lambda * trace * trace;
		energy = restVolume * psi;
	}

	// WARNING j = J* in c++ impl... used as a -> to array in code.
	static void computeGradCGreen(final float restVolume, final Matrix3 invRestMat, final Matrix3 sigma,
	                              Vector3[] j) {
		Matrix3 t = new Matrix3(invRestMat).inv();
		Matrix3 h = new Matrix3(sigma).mul(t).scl(restVolume);

		for (int i = 0; i < 3; i++) {
			j[i] = new Vector3(h.val[mtx3Pos(0, i)], h.val[mtx3Pos(1, i)], h.val[mtx3Pos(2, i)]);
		}
		j[3] = new Vector3(j[0]).sub(j[1]).sub(j[2]).scl(-1);
	}

	// note this is VERY similar to the computeGreenStrainAndPiolaStress method...
	static void computeGreenStrainAndPiolaStressInversion(
			final Vector3 x1, final Vector3 x2, final Vector3 x3, final Vector3 x4,
	        final Matrix3 invRestMat,
	        final float restVolume,
	        final float mu, final float lambda,
	        Matrix3 epsilon, Matrix3 sigma, float energy) {
		// Determine partial x/ partial m_1
		Matrix3 f = new Matrix3();
		final Vector3 p14 = new Vector3(x1).sub(x4);
		final Vector3 p24 = new Vector3(x2).sub(x4);
		final Vector3 p34 = new Vector3(x3).sub(x4);

		float[] irm = invRestMat.val;
		f.val[Matrix3.M00] = p14.x * irm[Matrix3.M00] + p24.x * irm[Matrix3.M10] + p34.x * irm[Matrix4.M20];
		f.val[Matrix3.M01] = p14.x * irm[Matrix3.M01] + p24.x * irm[Matrix3.M11] + p34.x * irm[Matrix4.M21];
		f.val[Matrix3.M02] = p14.x * irm[Matrix3.M02] + p24.x * irm[Matrix3.M12] + p34.x * irm[Matrix4.M22];
		f.val[Matrix3.M10] = p14.y * irm[Matrix3.M00] + p24.y * irm[Matrix3.M10] + p34.y * irm[Matrix4.M20];
		f.val[Matrix3.M11] = p14.y * irm[Matrix3.M01] + p24.y * irm[Matrix3.M11] + p34.y * irm[Matrix4.M21];
		f.val[Matrix3.M12] = p14.y * irm[Matrix3.M02] + p24.y * irm[Matrix3.M12] + p34.y * irm[Matrix4.M22];
		f.val[Matrix3.M20] = p14.z * irm[Matrix3.M00] + p24.z * irm[Matrix3.M10] + p34.z * irm[Matrix4.M20];
		f.val[Matrix3.M21] = p14.z * irm[Matrix3.M01] + p24.z * irm[Matrix3.M11] + p34.z * irm[Matrix4.M21];
		f.val[Matrix3.M22] = p14.z * irm[Matrix3.M02] + p24.z * irm[Matrix3.M12] + p34.z * irm[Matrix4.M22];

		Matrix3 ft_f = new Matrix3(f).transpose().mul(f);

		// inversion handling
		Matrix3 v = new Matrix3();
		Vector3 s = new Vector3();
		eigenDecomposition(ft_f, v, s); // fix references

		if (s.x < 0.0f) s.x = 0.0f;
		if (s.y < 0.0f) s.y = 0.0f;
		if (s.z < 0.0f) s.z = 0.0f;

		final float detV = v.det();
		if (detV < 0.0f) {
			float minLambda = Float.MAX_VALUE;
			int pos = 0;
			if (s.x < minLambda) { minLambda = s.x; }
			if (s.y < minLambda) { pos = 1; minLambda = s.y; }
			if (s.z < minLambda) { pos = 2; minLambda = s.z; }

			v.val[mtx3Pos(0, pos)] = -v.val[mtx3Pos(0, pos)];
			v.val[mtx3Pos(1, pos)] = -v.val[mtx3Pos(1, pos)];
			v.val[mtx3Pos(2, pos)] = -v.val[mtx3Pos(2, pos)];
		}

		Vector3 hatF = new Vector3((float)Math.sqrt(s.x), (float)Math.sqrt(s.y), (float)Math.sqrt(s.z));
		Matrix3 vt = new Matrix3(v).transpose();

		int chk = 0;
		int pos = 0;
		if (Math.abs(hatF.x) < 1.0e-4f) { chk++; }
		if (Math.abs(hatF.y) < 1.0e-4f) { pos = 1; chk++; }
		if (Math.abs(hatF.z) < 1.0e-4f) { pos = 2; chk++; }

		Matrix3 u;
		if (chk > 0) {
			if (chk > 1) u = new Matrix3();
			else {
				u = new Matrix3(f).mul(v);
				for (int i = 0; i < 3; i++) {
					if (i != pos) {
						for (int j = 0; j < 3; j++)
						u.val[mtx3Pos(j, i)] *= 1.0f / get(hatF, i);
					}
				}

				Vector3[] va = new Vector3[2];
				int index = 0;
				for (int i = 0; i < 3; i++) {
					if (i != pos) {
						va[index++] = new Vector3(u.val[mtx3Pos(0, i)], u.val[mtx3Pos(1, i)], u.val[mtx3Pos(2, i)]);
					}
				}

				Vector3 vec = va[0].crs(va[1]);
				vec.nor();
				u.val[mtx3Pos(0, pos)] = vec.x;
				u.val[mtx3Pos(1, pos)] = vec.y;
				u.val[mtx3Pos(2, pos)] = vec.z;
			}
		} else {
			Vector3 hatFInv = new Vector3(1.0f / hatF.x, 1.0f / hatF.y, 1.0f / hatF.z);
			u = new Matrix3(f).mul(v);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					u.val[mtx3Pos(j, i)] *= get(hatFInv, i);
				}
			}
		}

		final float detU = u.det();

		if (detU < 0.0f) {
			float minLambda = Float.MAX_VALUE;
			int _pos = 0;
			for (int i = 0; i < 3; i++) {
				if (hatF.x < minLambda) { minLambda = hatF.x; }
				if (hatF.y < minLambda) { _pos = 1; minLambda = hatF.y; }
				if (hatF.z < minLambda) { _pos = 2; minLambda = hatF.z; }
			}
			set(hatF, _pos, -get(hatF, _pos));
			u.val[mtx3Pos(0, _pos)] = -u.val[mtx3Pos(0, _pos)];
			u.val[mtx3Pos(1, _pos)] = -u.val[mtx3Pos(1, _pos)];
			u.val[mtx3Pos(2, _pos)] = -u.val[mtx3Pos(2, _pos)];
		}

		final float minXVal = 0.577f; //???
		for (int i = 0; i < 3; i++) {
			if (get(hatF, i) < minXVal) {
				set(hatF, i, minXVal);
			}
		}

		Vector3 epsilonHatF = new Vector3(
				0.5f * (hatF.x * hatF.x - 1.0f),
				0.5f * (hatF.y * hatF.y - 1.0f),
				0.5f * (hatF.z * hatF.z - 1.0f));

		final float trace = epsilonHatF.x + epsilonHatF.x + epsilonHatF.z;
		final float ltrace = lambda * trace;
		Vector3 sigmaVec = new Vector3(epsilonHatF).scl(2.0f * mu);
		sigmaVec.add(ltrace);
		sigmaVec.x = hatF.x * sigmaVec.x;
		sigmaVec.y = hatF.y * sigmaVec.y;
		sigmaVec.z = hatF.z * sigmaVec.z;

		Matrix3 sigmaDiag = new Matrix3();
		setRow(new Vector3(sigmaVec.x, 0, 0), sigmaDiag, 0);
		setRow(new Vector3(0, sigmaVec.y, 0), sigmaDiag, 1);
		setRow(new Vector3(0, 0, sigmaVec.z), sigmaDiag, 2);

		Matrix3 epsDiag = new Matrix3();
		setRow(new Vector3(epsilonHatF.x, 0, 0), epsDiag, 0);
		setRow(new Vector3(0, epsilonHatF.y, 0), epsDiag, 1);
		setRow(new Vector3(0, 0, epsilonHatF.z), epsDiag, 2);

		epsilon = new Matrix3(u).mul(vt);
		sigma = new Matrix3(u).mul(sigmaDiag).mul(vt);

		float psi = 0.0f;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				psi += epsilon.val[mtx3Pos(i, j)] * epsilon.val[mtx3Pos(i, j)];
			}
		}
		psi = mu * psi + 0.5f * lambda * trace * trace;
		energy = restVolume * psi;
	}

	// this looks like it could be optimized a bit.
	boolean solveFEMTetraConstraint(final Vector3 p0, final float invMass0,
	                                final Vector3 p1, final float invMass1,
	                                final Vector3 p2, final float invMass2,
	                                final Vector3 p3, final float invMass3,
	                                final float restVolume,
	                                final Matrix3 invRestMat,
	                                final float youngsModulus,
	                                final float poissonRatio,
	                                final boolean handleInversion,
	                                Vector3 corr0, Vector3 corr1, Vector3 corr2, Vector3 corr3) {
		if (youngsModulus <= 0.0f) {
			return true;
		}

		if (poissonRatio < 0.0f || poissonRatio > 0.49f) {
			return false;
		}

		corr0 = Vector3.Zero;
		corr1 = Vector3.Zero;
		corr2 = Vector3.Zero;
		corr3 = Vector3.Zero;
		float C = 0.0f;
		Vector3[] gradC = new Vector3[4];
		Matrix3 epsilon = new Matrix3();
		Matrix3 sigma = new Matrix3();
		float volume = new Vector3(p1).sub(p0).crs(new Vector3(p2).sub(p0)).dot(new Vector3(p3).sub(p0)) / 6.0f;
		float mu = youngsModulus / 2.0f / (1.0f + poissonRatio);
		float lambda = youngsModulus * poissonRatio / (1.0f + poissonRatio) / (1.0f - 2.0f * poissonRatio);

		if (!handleInversion || volume > 0.0f) {
			computeGreenStrainAndPiolaStress(p0, p1, p2, p3, invRestMat, restVolume, mu, lambda, epsilon, sigma, C);
		} else {
			computeGreenStrainAndPiolaStressInversion(p0, p1, p2, p3, invRestMat, restVolume, mu, lambda, epsilon, sigma, C);
		}
		computeGradCGreen(restVolume, invRestMat, sigma, gradC);

		float sumNormGradC = invMass0 * gradC[0].len2()
				+ invMass1 * gradC[1].len2()
				+ invMass2 * gradC[2].len2()
				+ invMass3 * gradC[3].len2();

		if (sumNormGradC < 1.0e-9f) {
			return false;
		}
		final float s = C / sumNormGradC;
		corr0 = gradC[0].scl(-s * invMass0);
		corr1 = gradC[1].scl(-s * invMass1);
		corr2 = gradC[2].scl(-s * invMass2);
		corr3 = gradC[3].scl(-s * invMass3);

		return true;
	}

	boolean computePBFDensity(final int particleIndex,
	                          final int numberOfParticles,
	                          final Vector3[] x,
	                          final float[] mass,
	                          final Vector3[] boundaryX,
	                          final float[] boundaryPsi,
	                          final int numNeighbors,
	                          final int[] neighbors,
	                          final float density0,
	                          final boolean boundaryHandling,
	                          float density_err, // pass by ref Float
	                          float density) {   // pass by ref Float
		density = mass[particleIndex] * sphKernel.wZero();
		for (int neighborIndex : neighbors) {
			if (neighborIndex < numberOfParticles) {
				density += mass[neighborIndex] * sphKernel.w(new Vector3(x[particleIndex]).sub(x[neighborIndex]));
			} else if (boundaryHandling) {
				density += boundaryPsi[neighborIndex - numberOfParticles] * sphKernel.w(new Vector3(x[particleIndex]).sub(boundaryX[neighborIndex - numberOfParticles]));
			}
		}
		density_err = Math.max(density, density0) - density0;
		return true;
	}

	boolean computePBFLagrangeMultiplier(final int particleIndex,
	                                     final int numberOfParticles,
	                                     final Vector3[] x,
	                                     final float[] mass,
	                                     final Vector3[] boundaryX,
	                                     final float[] boundaryPsi,
	                                     final float density,
	                                     final int numNeighbors,
	                                     final int[] neighbors,
	                                     final float density0,
	                                     final boolean boundaryHandling,
	                                     float lambda) {
		final float eps = 1.0e-6f;
		final float C = Math.max(density / density0 - 1.0f, 0f);

		if (C != 0f) {
			float sum_grad_C2 = 0f;
			Vector3 gradC_i = Vector3.Zero;

			for (int neighborIndex : neighbors) {
				if (neighborIndex < numberOfParticles) {
					Vector3 gradC_j = sphKernel.gradW(new Vector3(x[particleIndex]).sub(x[neighborIndex])).scl(-mass[neighborIndex] / density0);
					sum_grad_C2 += gradC_j.len2();
					gradC_i.sub(gradC_j);
				} else if (boundaryHandling) {
					Vector3 gradC_j = sphKernel.gradW(new Vector3(x[particleIndex]).sub(boundaryX[neighborIndex - numberOfParticles])).scl(-boundaryPsi[neighborIndex - numberOfParticles] / density0);
					sum_grad_C2 += gradC_j.len2();
					gradC_i.sub(gradC_j);
				}
			}
			sum_grad_C2 += gradC_i.len2();
			lambda = -C / (sum_grad_C2 + eps);
		}

		return true;
	}

	boolean solveDensityConstraint(final int particleIndex,
	                               final int numberOfParticles,
	                               final Vector3[] x,
	                               final float[] mass,
	                               final Vector3[] boundaryX,
	                               final float[] boundaryPsi,
	                               final int numNeighbors,
	                               final int[] neighbors,
	                               final float density0,
	                               final boolean boundaryHandling,
	                               final float[] lambda,
	                               Vector3 corr) {
		corr.setZero();
		for (int neighborIndex : neighbors) {
			if (neighborIndex < numberOfParticles) {
				Vector3 gradC_j = sphKernel.gradW(x[particleIndex].sub(x[neighborIndex])).scl(-mass[neighborIndex] / density0);
				corr.sub(gradC_j.scl(lambda[particleIndex] + lambda[neighborIndex]));
			} else if (boundaryHandling) {
				Vector3 gradC_j = sphKernel.gradW(x[particleIndex].sub(boundaryX[neighborIndex - numberOfParticles])).scl(-boundaryPsi[neighborIndex - numberOfParticles] / density0);
				corr.sub(gradC_j.scl(lambda[particleIndex]));
			}
		}
		return true;
	}
}
