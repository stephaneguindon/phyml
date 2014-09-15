package nz.org.nesi.phyml.swing;

import grisu.control.ServiceInterface;
import grisu.frontend.control.login.LoginException;
import grisu.frontend.control.login.LoginManager;
import grisu.frontend.view.swing.jobmonitoring.single.SingleJobTabbedPane;
import grisu.frontend.view.swing.login.GrisuLoginPanel;
import grisu.frontend.view.swing.login.LoginProgressPanel;
import grisu.frontend.view.swing.login.ServiceInterfaceHolder;
import grith.gridsession.GridClient;
import grith.jgrith.cred.AbstractCred;
import grith.jgrith.cred.Cred;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.swing.*;
import java.awt.*;

public class GridPanel extends JPanel implements ServiceInterfaceHolder {

	private static final long serialVersionUID = 8384764165781030979L;

	//public static final String APPLICATION_NAME = Constants.ALLJOBS_KEY;
    public static final String APPLICATION_NAME = "PhyML";

	private static final Logger myLogger = LoggerFactory.getLogger(GridPanel.class);

	private final String JOBLIST_PANEL = "ROOT";
	private final String LOGIN_PANEL = "LOGIN";
	private final String PROGRESS_PANEL = "PROGRESS";

	private ServiceInterface si;
	private GrisuLoginPanel grisuLoginPanel;
	private final GridClient gridClient;
	private LoginProgressPanel loginProgressPanel;
	private SingleJobTabbedPane simpleSingleJobsGrid;

	/**
	 * Create the panel.
	 */
	public GridPanel(GridClient gc) {
        gridClient = gc;

		setLayout(new CardLayout(0, 0));
		add(getGrisuLoginPanel(), LOGIN_PANEL);
		add(getProgressPanel(), PROGRESS_PANEL);

        new Thread() {
            public void run() {
                Cred c = AbstractCred.getExistingCredential();

                if (c != null && c.isValid()) {
                    ServiceInterface si = null;
                    try {
                        si = LoginManager.login("bestgrid", gridClient.getCredential(), false);
                        setServiceInterface(si);
                    } catch (LoginException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }

                }

            }
        }.start();


	}

	public ServiceInterface getServiceInterface() {
		return si;
	}

    public GridClient getGridClient() {
        return gridClient;
    }

	public void setServiceInterface(ServiceInterface si) {
		this.si = si;

		try {
			getProgressPanel().setLoginToBackend(si);
			switchToProgressPanel();
			add(getSimpleSingleJobsGrid(), JOBLIST_PANEL);

			switchToJobListPanel();
		} catch (final Exception ie) {
			myLogger.error(ie.getLocalizedMessage(), ie);
			switchToLoginPanel();
		}
	}

	private GrisuLoginPanel getGrisuLoginPanel() {
		if (grisuLoginPanel == null) {
			grisuLoginPanel = new GrisuLoginPanel(this);
			grisuLoginPanel.setSessionClient(gridClient);
		}
		return grisuLoginPanel;
	}
	private LoginProgressPanel getProgressPanel() {
		if (loginProgressPanel == null) {
			loginProgressPanel = new LoginProgressPanel();
		}
		return loginProgressPanel;
	}

	private void switchToProgressPanel() {

		SwingUtilities.invokeLater(new Thread() {
			@Override
			public void run() {
				final CardLayout cl = (CardLayout) (getLayout());
				cl.show(GridPanel.this, PROGRESS_PANEL);
			}
		});
	}

	private void switchToLoginPanel() {

		SwingUtilities.invokeLater(new Thread() {
			@Override
			public void run() {
				final CardLayout cl = (CardLayout) (getLayout());
				cl.show(GridPanel.this, LOGIN_PANEL);
			}
		});
	}


	private SingleJobTabbedPane getSimpleSingleJobsGrid() {
		if (simpleSingleJobsGrid == null) {
			simpleSingleJobsGrid = new SingleJobTabbedPane(si, APPLICATION_NAME);
		}
		return simpleSingleJobsGrid;
	}

	private void switchToJobListPanel() {

		SwingUtilities.invokeLater(new Thread() {
			@Override
			public void run() {
				final CardLayout cl = (CardLayout) (getLayout());
				cl.show(GridPanel.this, JOBLIST_PANEL);
			}
		});
	}
}
