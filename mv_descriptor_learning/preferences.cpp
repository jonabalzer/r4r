/*////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the R4R library.
//
// The R4R library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The R4R library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////*/

#include "preferences.h"
#include "ui_preferences.h"

#include <QFileDialog>

#include <iostream>


Preferences::Preferences(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Preferences),
    m_parent(parent) {

    ui->setupUi(this);

}

Preferences::~Preferences() {

    delete ui;

}

void Preferences::on_applyButton_clicked()
{

    // create new parameter object
    CParameters params;

    // write current parameter set to object
    params.Set("SCALE",ui->scaleSpinBox->value());
    params.Set("FEATURE_THRESHOLD",ui->thresholdEdit->text().toDouble());
    params.Set("MIN_NO_FEATURES",ui->minFeatEdit->text().toInt());
    params.Set("MAX_NO_FEATURES",ui->maxFeatEdit->text().toInt());
    params.Set("BUFFER_LENGTH",ui->ringBufferEdit->text().toInt());
    params.Set("TRACKING_HSIZE",ui->trackSizeEdit->text().toInt());
    params.Set("DESCRIPTOR_HSIZE",ui->descSizeEdit->text().toInt());
    params.Set("GRAD_SMOOTH_SIGMA",ui->sigmaEdit->text().toDouble());
    params.Set("MINIMAL_FEATURE_HDISTANCE_INIT",ui->minInitDistEdit->text().toInt());
    params.Set("MINIMAL_FEATURE_HDISTANCE_CLEAN",ui->minDistEdit->text().toInt());
    params.Set("MAX_HAMMING_DISTANCE",ui->hammingEdit->text().toInt());
    params.Set("LK_PYRAMID_LEVEL",ui->lkLevelSpinBox->value());
    params.Set("MAX_ITER",ui->maxIterEdit->text().toInt());
    params.Set("ACCURACY",ui->accEdit->text().toDouble());
    params.Set("LAMBDA",ui->lambdaEdit->text().toDouble());
    params.Set("ALPHA_GRAD_NORM",ui->alphaEdit->text().toDouble());
    params.Set("COMPUTE_ID",(int)ui->idOnCheckBox->isChecked());
    params.Set("NORMALIZE_ID",ui->idNormComboBox->currentIndex());
    params.Set("COMPUTE_GRAD",(int)ui->gradOnCheckBox->isChecked());
    params.Set("NORMALIZE_GRAD",ui->gradNormComboBox->currentIndex());
    params.Set("AGGREGATOR",ui->aggComboBox->currentIndex());
    params.Set("AGG_DS",ui->downSampleEdit->text().toInt());
    params.Set("COMPUTE_HOG",(int)ui->hogOnCheckBox->isChecked());

    // send this to main window
    emit params_changed(params);

    this->close();

}

void Preferences::on_saveButton_clicked() {

    // create new parameter object
    CParameters params;

    // write current parameter set to object
    params.Set("SCALE",ui->scaleSpinBox->value());
    params.Set("FEATURE_THRESHOLD",ui->thresholdEdit->text().toDouble());
    params.Set("MIN_NO_FEATURES",ui->minFeatEdit->text().toInt());
    params.Set("MAX_NO_FEATURES",ui->maxFeatEdit->text().toInt());
    params.Set("BUFFER_LENGTH",ui->ringBufferEdit->text().toInt());
    params.Set("TRACKING_HSIZE",ui->trackSizeEdit->text().toInt());
    params.Set("DESCRIPTOR_HSIZE",ui->descSizeEdit->text().toInt());
    params.Set("GRAD_SMOOTH_SIGMA",ui->sigmaEdit->text().toDouble());
    params.Set("MINIMAL_FEATURE_HDISTANCE_INIT",ui->minInitDistEdit->text().toInt());
    params.Set("MINIMAL_FEATURE_HDISTANCE_CLEAN",ui->minDistEdit->text().toInt());
    params.Set("MAX_HAMMING_DISTANCE",ui->hammingEdit->text().toInt());
    params.Set("LK_PYRAMID_LEVEL",ui->lkLevelSpinBox->value());
    params.Set("MAX_ITER",ui->maxIterEdit->text().toInt());
    params.Set("ACCURACY",ui->accEdit->text().toDouble());
    params.Set("LAMBDA",ui->lambdaEdit->text().toDouble());
    params.Set("ALPHA_GRAD_NORM",ui->alphaEdit->text().toDouble());
    params.Set("COMPUTE_ID",(int)ui->idOnCheckBox->isChecked());
    params.Set("NORMALIZE_ID",ui->idNormComboBox->currentIndex());
    params.Set("COMPUTE_GRAD",(int)ui->gradOnCheckBox->isChecked());
    params.Set("NORMALIZE_GRAD",ui->gradNormComboBox->currentIndex());
    params.Set("AGGREGATOR",ui->aggComboBox->currentIndex());
    params.Set("AGG_DS",ui->downSampleEdit->text().toInt());
    params.Set("COMPUTE_HOG",(int)ui->hogOnCheckBox->isChecked());

    QString filename = QFileDialog::getSaveFileName(this, tr("Save file..."),
                               ".",
                               tr("(*.txt)"));

    params.SaveToFile(filename.toStdString().c_str());


}

void Preferences::on_loadButton_clicked() {

    QString filename = QFileDialog::getOpenFileName(this, tr("Save file..."),
                                                    ".",
                                                    tr("(*.txt)"));

    CParameters params;
    params.OpenFromFile(filename.toStdString().c_str());

    ui->scaleSpinBox->setValue(params.GetIntParameter("SCALE"));
    ui->thresholdEdit->setText(QString::number(params.GetDoubleParameter("FEATURE_THRESHOLD")));
    ui->minFeatEdit->setText(QString::number(params.GetIntParameter("MIN_NO_FEATURES")));
    ui->maxFeatEdit->setText(QString::number(params.GetIntParameter("MAX_NO_FEATURES")));
    ui->ringBufferEdit->setText(QString::number(params.GetIntParameter("BUFFER_LENGTH")));
    ui->trackSizeEdit->setText(QString::number(params.GetIntParameter("TRACKING_HSIZE")));
    ui->descSizeEdit->setText(QString::number(params.GetIntParameter("DESCRIPTOR_HSIZE")));
    ui->sigmaEdit->setText(QString::number(params.GetDoubleParameter("GRAD_SMOOTH_SIGMA")));
    ui->minInitDistEdit->setText(QString::number(params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE_INIT")));
    ui->minDistEdit->setText(QString::number(params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE_CLEAN")));
    ui->hammingEdit->setText(QString::number(params.GetIntParameter("MAX_HAMMING_DISTANCE")));
    ui->lkLevelSpinBox->setValue(params.GetIntParameter("LK_PYRAMID_LEVEL"));
    ui->maxIterEdit->setText(QString::number(params.GetIntParameter("MAX_ITER")));
    ui->accEdit->setText(QString::number(params.GetDoubleParameter("ACCURACY")));
    ui->lambdaEdit->setText(QString::number(params.GetDoubleParameter("LAMBDA")));
    ui->alphaEdit->setText(QString::number(params.GetDoubleParameter("ALPHA_GRAD_NORM")));
    ui->idOnCheckBox->setChecked((bool)params.GetIntParameter("COMPUTE_ID"));
    ui->idNormComboBox->setCurrentIndex(params.GetIntParameter("NORMALIZE_ID"));
    ui->gradOnCheckBox->setChecked((bool)params.GetIntParameter("COMPUTE_GRAD"));
    ui->gradNormComboBox->setCurrentIndex(params.GetIntParameter("NORMALIZE_GRAD"));
    ui->aggComboBox->setCurrentIndex(params.GetIntParameter("AGGREGATOR"));
    ui->downSampleEdit->setText(QString::number(params.GetIntParameter("AGG_DS")));
    ui->hogOnCheckBox->setChecked((bool)params.GetIntParameter("COMPUTE_HOG"));

}

void Preferences::on_aggComboBox_currentIndexChanged(int index) {

    if(index==2)
        ui->downSampleEdit->setEnabled(true);
    else
        ui->downSampleEdit->setEnabled(false);

}
