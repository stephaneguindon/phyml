package phyml;

import java.awt.Color;

import javax.swing.JPanel;
/**
 * Implement all components for setting the "Branch support" variables.
 * 
 * @author Christoph Knapp
 *
 */
public class BranchSupport extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private BooStrAndAprLikRatTes booStrAndAprLikRatTes;

	/**
	 * Constructor to instantiate all components for setting the "Branch support" variables 
	 * and setting their size and location.
	 */
	public BranchSupport() {
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.25);
		add(new Separator("Branch Support", Color.lightGray, null,false));
		layout.setDimensions(1, 0.75);
		booStrAndAprLikRatTes = new BooStrAndAprLikRatTes();
		add(booStrAndAprLikRatTes);
	}

	/**
	 * Retrieves the number of bootstrapping trees from a booStrAndAprLikRatTes
	 * object.
	 * 
	 * @return String : number of bootstrapping trees converted into a String.
	 */
	public String getNumBootstrap() {
		return booStrAndAprLikRatTes.getNumBootstrap();
	}
}
