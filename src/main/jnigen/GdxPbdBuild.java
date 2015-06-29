import com.badlogic.gdx.jnigen.AntScriptGenerator;
import com.badlogic.gdx.jnigen.BuildConfig;
import com.badlogic.gdx.jnigen.BuildTarget;
import com.badlogic.gdx.jnigen.NativeCodeGenerator;

import java.io.File;

/**
 * The GdxPdbBuild generates the contents of the jni folder using the gdx jnigen tool.
 *
 * @since 1.1.0
 */
public class GdxPbdBuild {

	// Generate the files to the jni folder
	public static void main(String... args) throws Exception {
		new NativeCodeGenerator().generate("src/main/java", "target/classes", "jni");

		new AntScriptGenerator().generate(new BuildConfig("gdx-pbd"),
//				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.Windows, false), // win32
				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.Windows, true)   // win64
//				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.Linux, false),   // lin32
//				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.Linux, true),    // lin64
//				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.Android, false), // android
//				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.MacOsX, false),  // mac32
//				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.MacOsX, true),   // mac64
//				BuildTarget.newDefaultTarget(BuildTarget.TargetOs.IOS, false)      // ios
		);
	}
}
