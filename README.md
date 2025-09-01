# NMR 19F Spectrum Predictor
This application predicts the 19F NMR spectrum, the 19F NMR shifts for each fluorine atom in a molecule, and generates a molecular structure image based on a SMILES input string.


# Launching the Application
## 1. Start the Flask web application by running the following in the terminal:
    ```bash
    python app.py
    ```

You'll see output similar to this:
```bash
* Serving Flask app 'app'
* Debug mode: on
WARNING: This is a development server. Do not use it in a production deployment. Use a production WSGI server instead.
* Running on all addresses (0.0.0.0)
* Running on http://127.0.0.1:5002
* Running on http://192.168.1.180:5002
Press CTRL+C to quit
```

## 2. Copy the URL (e.g., http://192.168.1.180:5002) and open it in your web browser.
Using the Web Interface
On the webpage, input a SMILES string (a textual representation of chemical structure).
The application will display:
The predicted 19F NMR spectrum.
The 19F NMR shifts for each fluorine atom.
A molecular structure image.

### Notes
This application uses a development server and is not meant for production use. For production, use a WSGI server such as gunicorn or uwsgi.
If you encounter issues, ensure the required versions of Java and MATLAB are installed and accessible in your system's PATH.