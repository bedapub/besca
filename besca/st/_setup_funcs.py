from IPython.display import HTML
from IPython.display import Javascript


def create_button():
    input_form = """
    <div style="background-color:gainsboro; width:600px; padding:20px;">
    Continue running notebook with the given filtering parameters. <br>
    <button onclick="execute()">Continue</button>
    </div>
    """

    javascript = """
    <script type="text/Javascript">
        function execute(){
            var command = IPython.notebook.execute_cells_below();
            var kernel = IPython.notebook.kernel;
            kernel.execute(command);
        }
    """

    return HTML(input_form + javascript)


def create_popup():
    return Javascript(
        """
    function execute(){
            var command = IPython.notebook.execute_cells_below();
            var kernel = IPython.notebook.kernel;
            kernel.execute(command);
        }
    
    function interrupt(){
            var command = IPython.notebook.interrupt;
            var kernel = IPython.notebook.kernel;
            kernel.execute(command);
        }
    if(confirm("Continue running with the chosen standard parameters?")){execute()} else{interrupt()}
    """
    )
