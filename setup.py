from setuptools import setup, find_packages

setup(
    name="openmm_workflow",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
#    entry_points={
#        "console_scripts": [
#            "your_script_name = path.to.your.script:function_name",
#        ],
#    },
)

