package nz.org.nesi.phyml.swing;


import org.apache.commons.io.IOUtils;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;

public class PackageFileHelper {

	public static File TEMP_DIR = new File(
			System.getProperty("java.io.tmpdir"), "grisu-temp-files");

	public static File createTempFile(String content, String filename) {

		TEMP_DIR.mkdirs();

		File file = new File(TEMP_DIR, filename);
		file.deleteOnExit();

		try {
			InputStream is = new ByteArrayInputStream(content.getBytes());

			FileOutputStream fop = new FileOutputStream(file);

			LineEnds.convert(is, fop, LineEnds.STYLE_UNIX);
			fop.flush();
			fop.close();

			return file;
		} catch (Exception e3) {
			throw new RuntimeException(e3);
		}

	}

	public static File getFile(String fileName) {

		File file = new File(TEMP_DIR, fileName);
		file.deleteOnExit();

		if (!file.exists()) {
			file.getParentFile().mkdirs();
			InputStream inS = PackageFileHelper.class.getResourceAsStream("/"
					+ fileName);
			try {
				IOUtils.copy(inS, new FileOutputStream(file));
				inS.close();
			} catch (Exception e) {
				throw new RuntimeException(
						"Can't create temporary input files: "
								+ e.getLocalizedMessage());
			}
		}

		return file;
	}

	public static String getPath(String fileName) {

		return getFile(fileName).getAbsolutePath();
	}

}
