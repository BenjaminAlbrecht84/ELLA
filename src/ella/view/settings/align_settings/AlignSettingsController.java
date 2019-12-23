package ella.view.settings.align_settings;

import javafx.fxml.FXML;
import javafx.scene.control.*;

public class AlignSettingsController {

    @FXML
    protected Button runButton, cancelButton;

    @FXML
    protected Label info;

    // general-tab options ***************

    @FXML
    protected TextField reference, cores, multiplicity;

    @FXML
    protected TextArea inputArea, outputArea;

    @FXML
    protected ComboBox typeBox;

    @FXML
    protected Button refButton, inButton, outButton;

    // alignment-tab options *************

    @FXML
    protected ComboBox matrices, penalties;

    @FXML
    protected TextField frameshift;

    // output-tab options ****************

    @FXML
    protected Slider minCoverage;

    @FXML
    protected Label minCovLabel;

    @FXML
    protected TextField eValue;

}
