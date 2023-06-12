# pip install libsass watchdog
import sass
from http.server import SimpleHTTPRequestHandler, HTTPServer
from functools import partial
import mimetypes
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

def compile_scss():
    try:
        css = sass.compile(filename="./scss/custom.scss", output_style='compressed')
        with open("./wwwroot/custom.css", "w") as f:
            f.write(css)
        print("css updated")
    except Exception:
        pass

compile_scss()

class FileChangeEventHandler(FileSystemEventHandler):
    def on_any_event(self, event):
        compile_scss()

observer = Observer()
observer.schedule(FileChangeEventHandler(), "scss", recursive=False)
observer.start()

class RequestHandler(SimpleHTTPRequestHandler):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory="./wwwroot", **kwargs)

    if not mimetypes.inited:
        mimetypes.init() # try to read system mime.types
    extensions_map = mimetypes.types_map.copy()
    extensions_map.update({
        '': 'application/octet-stream', # Default
        '.py': 'text/plain',
        '.c': 'text/plain',
        '.h': 'text/plain',
        '.js': 'application/javascript',
        })    

handler_class = RequestHandler
httpd = HTTPServer(('', 8000), handler_class)

httpd.serve_forever()

observer.join()