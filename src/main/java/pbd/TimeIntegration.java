package pbd;

import com.badlogic.gdx.math.Vector3;

public class TimeIntegration {
	public static void semiImplicitEuler(final float h, final float mass,
	                                     Vector3 position, Vector3 velocity,
	                                     final Vector3 acceleration) {
		if (mass != 0.0f) {
			velocity.add(acceleration.scl(h));
			position.add(velocity.scl(h));
		}
	}

	public static void velocityUpdateFirstOrder(final float h, final float mass,
	                                            final Vector3 position, final Vector3 oldPosition,
	                                            Vector3 velocity) {
		if (mass != 0.0f) {
			velocity = new Vector3(position.sub(oldPosition).scl((1.0f / h)));
		}
	}

	public static void velocityUpdateSecondOrder(final float h, final float mass,
	                                             final Vector3 position, final Vector3 oldPosition,
	                                             final Vector3 positionOfLastStep, Vector3 velocity) {
		if (mass != 0.0f) {
			velocity = new Vector3(position.scl(1.5f).sub(oldPosition.scl(2f).add(positionOfLastStep.scl(0.5f))).scl(1.0f / h));
		}
	}
}
