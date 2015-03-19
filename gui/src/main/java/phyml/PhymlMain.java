package phyml;

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



//		Thread.setDefaultUncaughtExceptionHandler(new DefaultExceptionHandler());



        PhymlPanel frameContent = new PhymlPanel();
		@SuppressWarnings("unused")
		PhymlFrame startGui = new PhymlFrame("PhyML", 200, 100, 790, 820,
				JFrame.EXIT_ON_CLOSE, frameContent, true, true);
	}
}
