package phyml;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;

/**
 * Class CmdExec implements the call to phyml as a thread. This enables the
 * graphical user interface to get real time updates of the standard output.
 *
 * @author Christoph Knapp
 * @version 27-June-2012
 *
 */
public class CmdExec extends Thread {

	private String cmd;
	private String[] args;

	/**
	 * Constructor method for initialising an CmdExec object.
	 *
	 * @param cmd
	 *            A command for calling phyml from terminal.
	 */
	public CmdExec(String cmd) {
		this.cmd = cmd;
		args = new String[]{};
	}

	public CmdExec(String[] args) {
		this.cmd = "";
		this.args = args;
	}

	@Override
	public void run() {


		int exitStatus = -1;
		if (!cmd.equals("")&&args.length==0) {

			// trying to set executable flag, just in case.
			try {
			String token = cmd.split(" ")[0];

			File file = new File(token);
			file.setExecutable(true);
			} catch (Exception e) {
				e.printStackTrace();
			}

			try {
				Runtime rt = Runtime.getRuntime();
				Process process = rt.exec(cmd);
				BufferedReader input = new BufferedReader(
						new InputStreamReader(process.getInputStream()));
				String line1 = null;
				while ((line1 = input.readLine()) != null) {
					StandardOutPanel.setInput(line1);
					//Thread.sleep(20);
				}
				exitStatus = process.waitFor();
				System.out.println("Exit Status: " + exitStatus);
				PhymlPanel.SetSubmit(true);
				PhymlPanel.loadTrees();

			} catch (Throwable t) {
				t.printStackTrace();
			}
//		}else{
//			// basic housekeeping
//			LoginManager.initGrisuClient("phyml-grid");
//
//			// helps to parse commandline arguments, if you don't want to create
//			// your own parameter class, just use DefaultCliParameters
//			PhyMLParameters parameters = new PhyMLParameters();
//			// create the client
//			Client client = null;
//			try {
//				client = new Client(parameters, args);
//			} catch(Exception e) {
//				System.err.println("Could not start phyml-grid: "
//						+ e.getLocalizedMessage());
//				System.exit(1);
//			}
//
//			// finally:
//			// execute the "run"
//			client.run();
//			PhymlPanel.SetSubmit(true);
//			PhymlPanel.loadTrees();
		}
	}
}
