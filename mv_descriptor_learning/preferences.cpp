#include "preferences.h"
#include "ui_preferences.h"

#include <iostream>
#include <QFileDialog>

Preferences::Preferences(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Preferences),
    m_parent(parent)
{
    ui->setupUi(this);


}

Preferences::~Preferences()
{
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

    // send this to main window
    emit params_changed(params);

    this->close();

}

void Preferences::on_saveButton_clicked()
{

    // create new parameter object
    CParameters params;

    // write current parameter set to object
    params.Set("SCALE",ui->scaleSpinBox->value());
    params.Set("FEATURE_THRESHOLD",ui->thresholdEdit->text().toDouble());
    params.Set("MIN_NO_FEATURES",ui->minFeatEdit->text().toInt());
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

    QString filename = QFileDialog::getSaveFileName(this, tr("Save file..."),
                               ".",
                               tr("(*.txt)"));

    params.SaveToFile(filename.toStdString().c_str());


}

void Preferences::on_loadButton_clicked()
{

    QString filename = QFileDialog::getOpenFileName(this, tr("Save file..."),
                                                    ".",
                                                    tr("(*.txt)"));

    CParameters params;
    params.OpenFromFile(filename.toStdString().c_str());

    ui->scaleSpinBox->setValue(params.GetIntParameter("SCALE"));



}
