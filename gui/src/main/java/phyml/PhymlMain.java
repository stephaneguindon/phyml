package phyml;

import grisu.frontend.control.login.LoginManager;
import grisu.frontend.view.swing.utils.DefaultExceptionHandler;
import grisu.jcommons.utils.EnvironmentVariableHelpers;
import grith.gridsession.GridClient;
import grith.jgrith.Environment;

import javax.swing.*;

/**
 * Standard class implementing the main method which is necessary for all java
 * programs.
 *
 * @author Christoph Knapp
 * @date 02/07/12
 */
public class PhymlMain {

	/**
	 * Main method instantiating the JFrame of type PhymlFrame and the JPanel of
	 * type PhymlPanel.
	 *
	 * @param args
	 *            String[] : contains possible command line arguments. Unused in
	 *            this case.
	 */
	public static void main(String[] args) {

		// housekeeping for grid stuff
		LoginManager.initGrisuClient("phyml-grid-swing");

		LoginManager.setClientVersion(grisu.jcommons.utils.Version
				.get("this-client"));

		EnvironmentVariableHelpers.loadEnvironmentVariablesToSystemProperties();

		Thread.setDefaultUncaughtExceptionHandler(new DefaultExceptionHandler());

		Environment.initEnvironment();

        GridClient gc;
        try {
            gc = new GridClient();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        PhymlPanel frameContent = new PhymlPanel(gc);
		@SuppressWarnings("unused")
		PhymlFrame startGui = new PhymlFrame("PhyML", 200, 100, 790, 820,
				JFrame.EXIT_ON_CLOSE, frameContent, true, true);
	}
}
