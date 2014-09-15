package phyml;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.security.CodeSource;
import java.security.ProtectionDomain;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

/**
 * Class for unpacking files stored in a jar file.
 * 
 * @author Christoph Knapp
 * 
 */
public class UnpackTestDataFile {
	/**
	 * Starts u the process of extracting the file from the jar file.
	 * 
	 * @param positionToPut
	 *            String : full path + filename to where the file is extracted
	 *            to
	 * 
	 * @return String : returns the position to put variable.
	 * 
	 * @throws URISyntaxException
	 * @throws ZipException
	 * @throws IOException
	 */
	public static String start(String positionToPut) throws URISyntaxException,
			ZipException, IOException {
		final URI uri;
		uri = getJarURI();
		getFile(uri, positionToPut);
		return positionToPut;
	}

	/**
	 * Retrieves the URI to the jar file.
	 * 
	 * @return URI : position of the jar file.
	 * @throws URISyntaxException
	 */
	private static URI getJarURI() throws URISyntaxException {
		final ProtectionDomain domain;
		final CodeSource source;
		final URL url;
		final URI uri;

		domain = UnpackTestDataFile.class.getProtectionDomain();
		source = domain.getCodeSource();
		url = source.getLocation();
		uri = url.toURI();
		return (uri);
	}

	/**
	 * Retrieves the file from the jar file.
	 * 
	 * @param where
	 *            URI : where the jar file is located
	 * @param fileName
	 *            String : the full path to where the file is extracted to
	 *            including the file name.
	 * @throws ZipException
	 * @throws IOException
	 */
	private static void getFile(final URI where, final String fileName)
			throws ZipException, IOException {
		final File location;
		location = new File(where);
		// not in a JAR, just return the path on disk
		if (location.isDirectory()) {
			return;
		} else {
			final ZipFile zipFile;

			zipFile = new ZipFile(location);

			try {
				extract(zipFile, fileName);
			} finally {
				zipFile.close();
			}
		}
	}

	/**
	 * Extracts the specified file from the jar file.
	 * 
	 * @param zipFile
	 *            ZipFile : the jar file to be extracted from.
	 * @param fileName
	 *            String : the full path to where the file is extracted to
	 * 
	 * @throws IOException
	 */
	private static void extract(final ZipFile zipFile, final String fileName)
			throws IOException {
		final File tempFile;
		final ZipEntry entry;
		final InputStream zipStream;
		OutputStream fileStream;
		tempFile = new File(fileName);
		if (tempFile.exists()) {
			return;
		} else {
			if (!(new File(tempFile.getParent())).exists()) {
				(new File(tempFile.getParent())).mkdirs();
			}
			tempFile.createNewFile();
		}

		entry = zipFile.getEntry(tempFile.getName());

		if (entry == null) {
			throw new FileNotFoundException("cannot find file: " + fileName
					+ " in archive: " + zipFile.getName());
		}

		zipStream = zipFile.getInputStream(entry);
		fileStream = null;

		try {
			final byte[] buf;
			int i;

			fileStream = new FileOutputStream(tempFile);
			buf = new byte[1024];
			i = 0;

			while ((i = zipStream.read(buf)) != -1) {
				fileStream.write(buf, 0, i);
			}
		} finally {
			close(zipStream);
			close(fileStream);
		}
	}

	/**
	 * Closes a stream.
	 * 
	 * @param stream
	 *            Closeable : stream to close
	 */
	private static void close(final Closeable stream) {
		if (stream != null) {
			try {
				stream.close();
			} catch (final IOException ex) {
				ex.printStackTrace();
			}
		}
	}
}
