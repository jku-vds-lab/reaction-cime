import os

import uvicorn

dir = os.path.dirname(os.path.realpath(__file__))

if __name__ == "__main__":
    uvicorn.run(
        os.path.basename(dir) + ".dev_app:app",
        host="0.0.0.0",
        port=9000,
        reload=True,
        reload_dirs=[dir],
    )
