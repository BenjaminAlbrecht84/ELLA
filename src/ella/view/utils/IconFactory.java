/*
 * IconFactory.java Copyright (C) 2019. University of Tuebingen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package ella.view.utils;

import javafx.beans.property.ReadOnlyDoubleProperty;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;

import java.io.File;

public class IconFactory {

    public static ImageView createImageView(String name, ReadOnlyDoubleProperty size) {
        ImageView imageView = new ImageView(new Image(IconFactory.class.getResourceAsStream(File.separator +"ella"+ File.separator+"resources"+ File.separator + "icons" + File.separator + name)));
        imageView.fitHeightProperty().bind(size.multiply(0.7));
        imageView.fitWidthProperty().bind(size.multiply(0.7));
        return imageView;
    }

    public static ImageView createImageView(String name, int size) {
        ImageView imageView = new ImageView(new Image(IconFactory.class.getResourceAsStream(File.separator +"ella"+ File.separator+"resources"+ File.separator + "icons" + File.separator + name)));
        imageView.setFitWidth(size);
        imageView.setFitHeight(size);
        return imageView;
    }

}
