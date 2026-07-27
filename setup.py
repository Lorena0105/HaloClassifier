from setuptools import setup, find_packages
from pathlib import Path

# Read README.md for PyPI / documentation
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="HaloClassifier",
    version="1.0.1",
    packages=find_packages(),
    install_requires=[
        "pandas>=2.3.3,<3.0",
        "joblib>=1.5.2,<2.0",
        "tqdm>=4.67.1",
        "scikit-learn>=1.7.2,<1.8",
    ],
    entry_points={
        "console_scripts": [
            # Three aliases for the main menu
            "haloClassifier=haloClassifier_pkg.cli:menu",
            "haloclassifier=haloClassifier_pkg.cli:menu",
            "HaloClassifier=haloClassifier_pkg.cli:menu",

            # Operational commands
            "halo-generate-features=haloClassifier_pkg.cli:run_generate",
            "halo-classify=haloClassifier_pkg.cli:run_classify"
        ]
    },
    include_package_data=True,
    package_data={
        "haloClassifier_pkg": [
            "models/*",
            "models/*/*",
            "models/*/*/*"
        ]
    },

    # Include README.md
    long_description=long_description,
    long_description_content_type="text/markdown",

    python_requires=">=3.10",
)
