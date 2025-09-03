# 19F NMR Spectrum Predictor

A web application that predicts 19F NMR spectra for fluorinated compounds using machine learning. The application accepts SMILES strings as input and outputs predicted NMR spectra, molecular structures, and detailed analysis.

![19F NMR Predictor Web Application](web_page.png)

## 🌐 Live Application

**Try the application online**: [https://wqdx9jslij.execute-api.us-east-2.amazonaws.com/prod/](https://wqdx9jslij.execute-api.us-east-2.amazonaws.com/prod/)

Simply enter a SMILES string to get AI-powered 19F NMR predictions instantly.

### 🌟 Features

- **SMILES Input**: Accepts chemical structure notation for easy compound input
- **AI-Powered Prediction**: Uses a trained neural network for accurate 19F NMR spectrum prediction
- **Molecular Visualization**: Generates and displays molecular structures with fluorine atom highlighting
- **Confidence Levels**: Provides L1-L6 confidence indicators for each predicted peak
- **Professional UI**: Modern, responsive web interface built with Bootstrap 5
- **Real-time Validation**: Client-side SMILES validation with immediate feedback

### 📊 Performance

- **Prediction Time**: 50ms to 1 seconds (depending on compound complexity)

### 🧪 How to Use

#### 1. **Enter SMILES String**
Input a SMILES representation of your fluorinated compound. Examples:
- **PFBA**: `C(=O)(C(C(C(F)(F)F)(F)F)(F)F)O`
- **PFBS**: `C(C(C(F)(F)S(=O)(=O)O)(F)F)(C(F)(F)F)(F)F`

#### 2. **Get Predictions**
The application will generate:
- **Molecular Structure**: Visual representation with fluorine atoms highlighted in gold
- **Predicted 19F NMR Spectrum**: Peak positions with confidence level labels (L1-L6)
- **Detailed Analysis**: Confidence levels and summary statistics

## 🔧 Technical Details

### Architecture
- **Backend**: Flask web framework
- **ML Model**: TensorFlow Sequential neural network
- **Chemistry**: RDKit for molecular handling and visualization
- **Visualization**: Matplotlib for NMR spectrum generation
- **Frontend**: Bootstrap 5, Font Awesome, custom CSS

### Model Information
- **Type**: Feed-forward neural network (FFNN)
- **Input**: Molecular descriptors
- **Output**: 19F NMR chemical shifts
- **Training**: Fluorinated compounds dataset
- **Performance**: Optimized for speed and accuracy

## 🚀 Local Deployment

### Option 1: Docker (Recommended)

```bash
# Pull the image from Docker Hub
docker pull dandanrao/19nmr-predictor_ffnn:latest

# Run the container
docker run -d -p 5002:5002 --name nmr-predictor dandanrao/19nmr-predictor_ffnn:latest

# Access the web app
open http://localhost:5002
```

### Option 2: Clone the repo

```bash
# Clone the repository
git clone https://github.com/Dandan-Rao/Application_19F_NMR_Spectra_predictor_FFNN.git
cd Application_19F_NMR_Spectra_predictor_FFNN

# Install dependencies
pip install -r requirements.txt

# Run the application
python app.py

# Access the web app
open http://localhost:5002
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- **RDKit**: Chemical informatics toolkit
- **TensorFlow**: Machine learning framework
- **Flask**: Web framework
- **Bootstrap**: UI framework

## 📞 Support

For questions or issues:
- Create an issue in the repository
- Check the troubleshooting section
- Review the logs for error details

---

**Note**: This application is designed for research and educational purposes. Always verify predictions with experimental data for critical applications.
