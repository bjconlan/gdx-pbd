package pbd;

import com.badlogic.gdx.math.Vector3;

public class TimeIntegration {
	public static void semiImplicitEuler(final float h, final float mass,
	                                     Vector3 position, Vector3 velocity,
	                                     final Vector3 acceleration) {}

	public static void velocityUpdateFirstOrder(final float h, final float mass,
	                                            final Vector3 position, final Vector3 oldPosition,
	                                            Vector3 velocity) {}

	public static void velocityUpdateSecondOrder(final float h, final float mass,
	                                             final Vector3 position, final Vector3 oldPosition,
	                                             final Vector3 positionOfLastStep, Vector3 velocity) {
	}
}
