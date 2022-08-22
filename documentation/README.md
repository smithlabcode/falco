# falco documentation

This is the documentation for falco that uses
[mkdocs](https://mkdocs.readthedocs.io) to generate readthedocs pages.
The public web verison of this documentation is available at
[falco.readthedocs.io](https://falco.readthedocs.io), but for
uses who wish to see the documentation on a web browser offline, you
can build the documentation locally as described below.

### Dependencies

To build the documentation locally, install mkdocs

```
pip install -U mkdocs
```

### Local compilation

Build the HTML documentation by running
```
mkdocs build
```
which will create a `site` directory where markdown files are
converted to HTML

Create a local host for the HTML documentation by running

```
mkdocs serve
```

This will create the documentation, usually at http://localhost:8000 .
