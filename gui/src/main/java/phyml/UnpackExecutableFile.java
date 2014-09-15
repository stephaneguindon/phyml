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
 * Class to unpack a file from a jar file and put it a specified location.
 * 
 * @author Christoph Knapp
 * 
 */
public class UnpackExecutableFile {
	/**
	 * Starts up the process of extracting the file.
	 * 
	 * @param positionToPut
	 *            String : full path to were the file is extracted to.
	 * @return URI : with the position were the file is extracted to
	 * @throws URISyntaxException
	 * @throws ZipException
	 * @throws IOException
	 */
	public static URI start(String positionToPut) throws URISyntaxException,
			ZipException, IOException {
		final URI uri;
		final URI exe;
		uri = getJarURI();
		exe = getFile(uri, positionToPut);
		return exe;
	}

	/**
	 * Retrieves the position of the jar file.
	 * 
	 * @return URI : contains the position of the jar file.
	 * @throws URISyntaxException
	 */
	private static URI getJarURI() throws URISyntaxException {
		final ProtectionDomain domain;
		final CodeSource source;
		final URL url;
		final URI uri;
		domain = UnpackExecutableFile.class.getProtectionDomain();
		source = domain.getCodeSource();
		url = source.getLocation();
		uri = url.toURI();
		return (uri);
	}

	/**
	 * Tests whether the file is in the archive and starts the process of
	 * extracting it.
	 * 
	 * @param where
	 *            URI : the position of the jar file.
	 * @param fileName
	 *            String : position to put the file.
	 * @return URI the position where the file was put.
	 * @throws ZipException
	 * @throws IOException
	 */
	private static URI getFile(final URI where, final String fileName)
			throws ZipException, IOException {
		final File location;
		final URI fileURI;

		location = new File(where);

		// not in a JAR, just return the path on disk
		if (location.isDirectory()) {
			fileURI = URI.create(fileName);
		} else {
			final ZipFile zipFile;

			zipFile = new ZipFile(location);

			try {
				fileURI = extract(zipFile, fileName);
			} finally {
				zipFile.close();
			}
		}

		return (fileURI);
	}

	/**
	 * Extracts the file from the jar file.
	 * 
	 * @param zipFile
	 *            ZipFile : position of the jar file.
	 * @param fileName
	 *            String : position where to put the file.
	 * @return URI : the position where the file was extracted to
	 * 
	 * @throws IOException
	 */
	private static URI extract(final ZipFile zipFile, final String fileName)
			throws IOException {
		final File tempFile;
		final ZipEntry entry;
		final InputStream zipStream;
		OutputStream fileStream;

		tempFile = new File(fileName);
		if (tempFile.exists()) {
			return (tempFile.toURI());
		} else {
			if (!(new File(tempFile.getParent())).exists()) {
				File f = new File(tempFile.getParent());
				f.mkdirs();
				f.deleteOnExit();
			}
			tempFile.createNewFile();
			tempFile.setExecutable(true, false);
		}
		tempFile.deleteOnExit();
		tempFile.getParentFile().deleteOnExit();

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

		return (tempFile.toURI());
	}

	/**
	 * Closes the zipStream.
	 * 
	 * @param stream
	 *            Closeable : a closeable stream to close.
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
