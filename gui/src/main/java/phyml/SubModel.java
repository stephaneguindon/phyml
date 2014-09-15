package phyml;

import javax.swing.*;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.File;

/**
 * Class for instantiating all components to specify the actual substitution
 * model as well as setting there size and location.
 *
 * @author Christoph Knapp
 */

public class SubModel extends JPanel implements ActionListener {
    /**
     * default id
     */
    private static final long serialVersionUID = 1L;
    private JComboBox modelBox;
    private JLabel sMLab;
    private JLabel cMLab;
    private JLabel eFLab;
    private JLabel oRPLab;
    private CustomTextField curTextField;
    private JComboBox equiBox;
    private FreqPanel freqPanel;
    private String modelPath;
    private String molecularType;
    private JRadioButton optimiseRateParameterYes;
    private JRadioButton optimiseRateParameterNo;
    private OptimiseFreqAndRatio optimiseFreqAndRatioPanel;
    private ProportionAndSubrate proportionAndSubratePanel;
    private NumSubCatGammaAverage numSubCatGammaAveragePanel;

    /**
     * Constructor method for instantiating all components for specifying the
     * substitution model and setting there size and location.
     *
     * @param molecularType Setting either "DNA" as molecule type or "AA".
     */
    public SubModel(String molecularType) {
        this.molecularType = molecularType;
        modelPath = "";
        curTextField = new CustomTextField();
        String[] models;
        String[] equiChoices;
        if (molecularType.equals("DNA")) {
            models = new String[]{"HKY85", "F84", "TN93", "GTR", "custom",
                    "JC69", "K80"};
            equiChoices = new String[]{"empirical", "optimised", "user defined"};
            eFLab = new JLabel("Equilibrium Frequency");
        } else {
            models = new String[]{"LG", "WAG", "Dayhoff", "JTT", "Blossum62",
                    "Mt Rev", "Rt Rev", "Cp Rev", "DcMut", "VT", "Mt Mam",
                    "Mt Art", "HIVw", "HIVb", "Read from file"};
            equiChoices = new String[]{"model", "empirical"};
            eFLab = new JLabel("Amino Acid Frequency");
        }
        modelBox = new JComboBox(new DefaultComboBoxModel(models));
        modelBox.addActionListener(this);
        equiBox = new JComboBox(equiChoices);
        equiBox.addActionListener(this);
        optimiseRateParameterYes = new JRadioButton("Yes");
        optimiseRateParameterNo = new JRadioButton("No");
        optimiseRateParameterYes.addActionListener(this);
        optimiseRateParameterNo.addActionListener(this);
        setOptimiseRateOff(true);
        setOptimiseRate("NO");
        freqPanel = new FreqPanel();
        freqPanel.setCompVisible(true);
        sMLab = new JLabel("Substitution Model");
        cMLab = new JLabel("Current Model");
        oRPLab = new JLabel("Optimise Rate Parameter");
        if (molecularType.equals("AA")) {
            setCompVisible(false);
        }

        CustomGridLayout layout = new CustomGridLayout();
        setLayout(layout);
        layout.setDimensions(1, 0.05);
        add(new JPanel());
        layout.setDimensions(0.01, 0.15);
        add(new JPanel());
        layout.setDimensions(0.33, 0.15);
        add(sMLab);
        layout.setDimensions(0.27, 0.15);
        add(modelBox);
        layout.setDimensions(0.39, 0.15);
        add(new JPanel());
        layout.setDimensions(1, 0.05);
        add(new JPanel());
        layout.setDimensions(0.01, 0.15);
        add(new JPanel());
        layout.setDimensions(0.33, 0.15);
        add(cMLab);
        layout.setDimensions(0.27, 0.15);
        add(curTextField);
        layout.setDimensions(0.39, 0.15);
        add(new JPanel());
        layout.setDimensions(1, 0.05);
        add(new JPanel());
        layout.setDimensions(0.01, 0.15);
        add(new JPanel());
        layout.setDimensions(0.33, 0.15);
        add(eFLab);
        layout.setDimensions(0.27, 0.15);
        add(equiBox);
        layout.setDimensions(0.39, 0.15);
        add(new JPanel());
        layout.setDimensions(1, 0.05);
        add(new JPanel());
        layout.setDimensions(0.01, 0.15);
        add(new JPanel());
        layout.setDimensions(0.33, 0.15);
        add(oRPLab);
        layout.setDimensions(0.27, 0.15);
        add(new JPanel());
        layout.setDimensions(0.38, 0.15);
        JPanel p1 = new JPanel();
        CustomGridLayout lo1 = new CustomGridLayout();
        p1.setLayout(lo1);
        add(p1);
        lo1.setDimensions(0.2, 1);
        p1.add(new JPanel());
        lo1.setDimensions(0.4, 1);
        p1.add(optimiseRateParameterYes);
        lo1.setDimensions(0.4, 1);
        p1.add(optimiseRateParameterNo);
        layout.setDimensions(1, 0.05);
        add(new JPanel());
        layout.setDimensions(1, 0.15);
        add(freqPanel);
        setSelectedCompEnabled(2);
    }

    public void setProportionAndSubratePanel(ProportionAndSubrate proportionAndSubratePanel) {
        this.proportionAndSubratePanel = proportionAndSubratePanel;
    }

    public void setOptimiseFreqAndRatioPanel(OptimiseFreqAndRatio optimiseFreqAndRatioPanel) {
        this.optimiseFreqAndRatioPanel = optimiseFreqAndRatioPanel;
    }

    /**
     * This method sets all components in the gui to visible or not visible.
     *
     * @param b boolean : true if visible, false if not visible.
     */
    private void setCompVisible(boolean b) {
        cMLab.setVisible(b);
        oRPLab.setVisible(b);
        curTextField.setVisible(b);
        optimiseRateParameterYes.setVisible(b);
        optimiseRateParameterNo.setVisible(b);
    }

    /**
     * Changes the Comboboxmodel for the substitution model dropdown menu to DNA
     * models or AA models.
     *
     * @param molType String : either "DNA" or "AA".
     */
    public void setMoleculeType(String molType) {
        molecularType = molType;
        if (molecularType.equals("DNA")) {
            modelBox.setModel(new DefaultComboBoxModel(new String[]{"HKY85",
                    "F84", "TN93", "GTR", "custom", "JC69", "K80"}));
            setCompVisible(true);
            equiBox.setModel(new DefaultComboBoxModel(new String[]{"empirical",
                    "optimised", "user defined"}));
            eFLab.setText("Equilibrium Frequency");
            freqPanel.setVisible(true);
        } else if (molecularType.equals("AA")) {
            modelBox.setModel(new DefaultComboBoxModel(new String[]{"LG",
                    "WAG", "Dayhoff", "JTT", "Blossum62", "Mt Rev", "Rt Rev",
                    "Cp Rev", "DcMut", "VT", "Mt Mam", "Mt Art", "HIVw",
                    "HIVb", "Read from file"}));
            setCompVisible(false);
            equiBox.setModel(new DefaultComboBoxModel(new String[]{"model",
                    "empirical"}));
            eFLab.setText("Amino Acid Frequency");
            freqPanel.setVisible(false);
        }
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        if (e.getSource() == modelBox) {
            if (modelBox.getSelectedItem().toString().equals("JC69") ||
                    modelBox.getSelectedItem().toString().equals("K80")) {
                setOptimiseRateOff(true);
                setOptimiseRate("NO");
                PhymlPanel.setRatioBoxOff(false);
                equiBox.setSelectedIndex(2);
                setFrequencies("0.25", "0.25", "0.25", "0.25");
                setSelectedCompEnabled(1);

            } else if (modelBox.getSelectedItem().toString().equals("custom")) {
                setOptimiseRateOff(true);
                setOptimiseRate("YES");
                PhymlPanel.setRatioBoxOff(false);
                setSelectedCompEnabled(true);
                setFrequencies("", "", "", "");
            } else if (modelBox.getSelectedItem().toString().equals("Read from file")) {
                setOptimiseRateOff(true);
                setOptimiseRate("NO");
                JFileChooser fc = new JFileChooser();
                fc.setCurrentDirectory(new File(System.getProperty("user.dir")));
                int returnVal = fc.showOpenDialog(this);
                if (returnVal == JFileChooser.APPROVE_OPTION) {
                    modelPath = fc.getSelectedFile().getAbsolutePath();
                }
                if (modelPath.equals("")) {
                    modelBox.setSelectedIndex(0);
                }
                setFrequencies("", "", "", "");
            } else if (modelBox.getSelectedItem().toString().equals("GTR")) {
                setOptimiseRateOff(false);
                setOptimiseRate("NO");
                PhymlPanel.setRatioBoxOff(true);
                PhymlPanel.setRatioBoxDefault();
            } else {
                setOptimiseRateOff(true);
                setOptimiseRate("NO");
                PhymlPanel.setRatioBoxOff(false);
                modelPath = "";
                equiBox.setSelectedIndex(0);
                setSelectedCompEnabled(2);
                setFrequencies("", "", "", "");
            }
        } else if (e.getSource() == equiBox) {
            if (equiBox.getSelectedItem().toString().equals("user defined")) {
                freqPanel.setCompEnabled(true);
                freqPanel.setCompVisible(true);
            } else {
                freqPanel.setCompEnabled(false);
                freqPanel.setCompVisible(true);
            }
        }
        if (e.getSource() == optimiseRateParameterYes) {
            if (optimiseRateParameterYes.isSelected()) {
                optimiseRateParameterNo.setSelected(false);
                optimiseFreqAndRatioPanel.setRatioBoxOff(false);
            } else {
                optimiseRateParameterNo.setSelected(true);
                optimiseFreqAndRatioPanel.setRatioBoxOff(true);
            }
            setOptimiseRate("YES");
        }
        if (e.getSource() == optimiseRateParameterNo) {
            if (optimiseRateParameterNo.isSelected()) {
                optimiseRateParameterYes.setSelected(false);
                optimiseFreqAndRatioPanel.setRatioBoxOff(true);
            } else {
                optimiseRateParameterYes.setSelected(true);
                optimiseFreqAndRatioPanel.setRatioBoxOff(false);
            }
            setOptimiseRate("NO");
        }
    }

    private void setOptimiseRate(String string) {
        if (string.equals("YES")) {
            optimiseRateParameterYes.setSelected(true);
            optimiseRateParameterNo.setSelected(false);

            if ( optimiseFreqAndRatioPanel != null ) {
                optimiseFreqAndRatioPanel.setRatioBoxOff(false);
            }
            if ( proportionAndSubratePanel != null ) {
                proportionAndSubratePanel.disableTransationTransversionRatio(false);
            }

            if (numSubCatGammaAveragePanel != null) {
                numSubCatGammaAveragePanel.disableGammaShape(false);
            }

        } else if (string.equals("NO")) {
            optimiseRateParameterYes.setSelected(false);
            optimiseRateParameterNo.setSelected(true);

            if ( optimiseFreqAndRatioPanel != null ) {
                optimiseFreqAndRatioPanel.setRatioBoxOff(true);
            }

            if ( proportionAndSubratePanel != null ) {
                proportionAndSubratePanel.setProportionInvariableSitesToFixed(true);
                proportionAndSubratePanel.disableTransationTransversionRatio(true);
            }

            if (numSubCatGammaAveragePanel != null) {
                numSubCatGammaAveragePanel.disableGammaShape(true);
                numSubCatGammaAveragePanel.setGammaShapeFixed(true);
            }

        }
    }

    private void setOptimiseRateOff(boolean b) {
        optimiseRateParameterYes.setEnabled(!b);
        optimiseRateParameterNo.setEnabled(!b);
    }

    /**
     * Enables or disables components which should only be enabled if a user
     * defined model is choosen.
     *
     * @param i int : which combination of enabled components to select.
     */
    private void setSelectedCompEnabled(int i) {
        if (i == 1) {
            curTextField.setEnabled(false);
            equiBox.setEnabled(false);
            freqPanel.setCompEnabled(false);
        } else if (i == 2) {
            curTextField.setEnabled(false);
            equiBox.setEnabled(true);
            freqPanel.setCompEnabled(false);
        }
    }

    /**
     * Enables or disables components which should only be enabled if a user
     * defined model is choosen.
     *
     * @param b boolean : true if component enabled, false if disabled.
     */
    private void setSelectedCompEnabled(boolean b) {
        curTextField.setEnabled(b);
        equiBox.setEnabled(b);
    }

    public void setFrequencies(String freqA, String freqC, String freqG, String freqT) {
        freqPanel.setFrequencies(freqA, freqC, freqG, freqT);
    }

    /**
     * The selected model of substitution.
     *
     * @return String : The selected mode i.e "HKY85". If custom is selected the
     *         String represantation of the model consisting of 6 digits is
     *         returned i.e. "000000".
     */
    public String getSubModel() {
        if (molecularType.equals("DNA")
                && modelBox.getSelectedItem().toString().equals("custom")) {
            return curTextField.getText();
        }
        return modelBox.getSelectedItem().toString();
    }

    /**
     * If the moleeculetype AA is selected the user has the possibility to
     * provide a rate file.
     *
     * @return String : The full path to a specified Amino acid rate file.
     */
    public String getAARateFileName() {
        return modelPath;
    }

    /**
     * Retrieves the command line options for the specified model.
     *
     * @return String : Comand line option of substitution model frequencies. -f
     *         e, m, or “fA,fC,fG,fT” Nucleotide or amino-acid frequencies. – e
     *         : the character frequencies are determined as follows : -
     *         Nucleotide sequences: (Empirical) the equilibrium base
     *         frequencies are estimated by counting the occurence of the
     *         different bases in the alignment. - Amino-acid sequences:
     *         (Empirical) the equilibrium amino-acid frequencies are estimated
     *         by counting the occurence of the different amino-acids in the
     *         alignment. – m : the character frequencies are determined as
     *         follows : - Nucleotide sequences: (ML) the equilibrium base
     *         frequencies are estimated using maximum likelihood. - Amino-acid
     *         sequences: (Model) the equilibrium amino-acid frequencies are
     *         estimated using the frequencies defined by the substitution
     *         model. – “fA,fC,fG,fT” : only valid for nucleotide-based models.
     *         fA, fC, fG and fT are floating numbers that correspond to the
     *         frequencies of A, C, G and T respectively.
     */
    public String getFrequencies() {
        if (molecularType.equals("DNA")) {
            if (equiBox.getSelectedItem().toString().equals("empirical")) {
                return "e";
            } else if (equiBox.getSelectedItem().toString().equals("optimised")) {
                return "m";
            } else {
                return formatFrequencies(freqPanel.getFreq());
            }
        } else {
            if (equiBox.getSelectedItem().toString().equals("model")) {
                return "m";
            } else {
                return "e";
            }
        }
    }

    /**
     * Formats the frequencies for the custom substitution model.
     *
     * @param freq String[] : array of frequencies.
     * @return String : The frequencies in form of a comma seperated list for
     *         passing them into the command line.
     */
    private String formatFrequencies(String[] freq) {
        return freq[0] + "," + freq[1] + "," + freq[2] + "," + freq[3];
    }

    /**
     * Whether the Rate parameter needs to be optimised.
     *
     * @return boolean : true if the rate parameter should be optimised, false
     *         otherwise.
     */
    public boolean isOptimisedRateParameter() {
        if (optimiseRateParameterYes.isSelected()) {
            return true;
        }
        return false;
    }

    public void setNumSubCatGammaAveragePanel(NumSubCatGammaAverage nSCGA) {
        numSubCatGammaAveragePanel = nSCGA;
    }

    class CustomTextField extends JTextField implements ActionListener,
            FocusListener {
        /**
         * default id
         */
        private static final long serialVersionUID = 1L;

        public CustomTextField() {
            super("000000", 6);
            this.addActionListener(this);
            this.addFocusListener(this);
            this.setMinimumSize(new Dimension(40, 4));
            this.setMaximumSize(new Dimension(80, 25));
            this.setPreferredSize(new Dimension(60, 18));
        }

        @Override
        protected Document createDefaultModel() {
            return new NumOnlyDocument();
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            if (this.getText().length() < 6) {
                String txt = this.getText();
                for (int i = this.getText().length(); i < 6; i++) {
                    txt = txt + "0";
                }
                this.setText(txt);
            }
        }

        @Override
        public void focusGained(FocusEvent e) {
        }

        @Override
        public void focusLost(FocusEvent e) {
            if (this.getText().length() < 6) {
                String txt = this.getText();
                for (int i = this.getText().length(); i < 6; i++) {
                    txt = txt + "0";
                }
                this.setText(txt);
            }
        }

        /**
         * @author Christoph Knapp
         */
        class NumOnlyDocument extends PlainDocument {
            /**
             * default id
             */
            private static final long serialVersionUID = 1L;

            @Override
            public void insertString(int offs, String str, AttributeSet a)
                    throws BadLocationException {
                if (str == null || offs > 5 || !isNumber(str)
                        || this.getLength() > 5) {
                    return;
                }
                char[] upper = str.toCharArray();
                super.insertString(offs, new String(upper), a);
            }

            private boolean isNumber(String str) {
                try {
                    Integer.parseInt(str);
                    return true;
                } catch (NumberFormatException e) {
                    return false;
                }
            }
        }
    }
}
