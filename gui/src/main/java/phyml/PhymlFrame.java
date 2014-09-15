package phyml;

import java.awt.Container;
import java.awt.Dimension;

import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * Highly customisable standard frame for graphical user interface.
 * 
 * @author Christoph Knapp
 * @date 29/06/12
 */
public class PhymlFrame extends JFrame {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Highly customisable standard JFrame for graphical user interface.
	 * 
	 * @param title
	 *            title of the gui displayed on the bar on top of the gui
	 * @param x
	 *            : The initial y position of the gui on the screen (upper left
	 *            corner)
	 * @param y
	 *            : The initial y position of the gui on the screen (upper left
	 *            corner)
	 * @param width
	 *            : The initial width of the JPanel frameContent
	 * @param height
	 *            : The initial height of the JPanel frameContent
	 * @param closOP
	 *            : Integer specifying what to do when the gui is closed i.e.
	 *            JFrame.EXIT_ON_CLOSE
	 * @param frameContent
	 *            : JPanel to be added to the JFrame
	 * @param vis
	 *            : boolean specifying wether the gui is visible or not
	 * @param resize
	 *            : booloean value specifying wether the gui is resizable or not
	 * 
	 */
	public PhymlFrame(String title, int x, int y, int width, int height,
			int closeOP, JPanel frameContent, boolean vis, boolean resize) {
		setTitle(title);
		setLocation(x, y);
		setDefaultCloseOperation(closeOP);
		this.setResizable(resize);
		Container visibleArea = getContentPane();
		visibleArea.add(frameContent);
		frameContent.setPreferredSize(new Dimension(width, height));
		pack();
		frameContent.requestFocusInWindow();
		setVisible(vis);
	}
}
