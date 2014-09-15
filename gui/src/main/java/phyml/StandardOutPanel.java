package phyml;

import grisu.jcommons.utils.swing.SmartScroll;

import javax.swing.*;
import javax.swing.text.DefaultCaret;
import java.awt.*;

/**
 * Extends JPanel and implements the components to display the standard output
 * of phyml.
 *
 * @author Christoph Knapp
 *
 */

public class StandardOutPanel extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private static JTextArea editorPane;
	private static  JScrollPane editorScrollPane;

	/**
	 * Constructor method implements the components to display the standard
	 * output of phyml.
	 */
	public StandardOutPanel() {
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 1);
		editorPane = new JTextArea();
		DefaultCaret caret = (DefaultCaret)editorPane.getCaret();
		caret.setUpdatePolicy(DefaultCaret.ALWAYS_UPDATE);
		Font font = new Font("Courier", Font.PLAIN, 12);
		editorPane.setFont(font);
		editorPane.setEditable(false);
		editorPane.setDropMode(DropMode.INSERT);
		editorScrollPane = new JScrollPane(editorPane);
        new SmartScroll(editorScrollPane);
		editorScrollPane
				.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		editorScrollPane
				.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		add(editorScrollPane);
	}

	/**
	 * Adds a new line to the standard output panel.
	 *
	 * @param line
	 *            String : A line to add to the output panel
	 */
	public static void setInput(String line) {
		if(line.contains("\t")){
			System.out.println(line);
		}
		editorPane.append(line + "\n");
	}


    public static void clearPanel() {

        editorPane.setText("");

    }
}
