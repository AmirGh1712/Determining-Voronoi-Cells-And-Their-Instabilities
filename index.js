const express = require('express')
const bodyParser = require('body-parser');
const path = require('path')
const PythonShell = require('python-shell').PythonShell;

// Setting the server.
const app = express()
const PORT = process.env.PORT || 3000
const server = app.listen(PORT, function() {
    console.log('listening on ' + PORT)
})
app.use(express.static('public'))
app.use(bodyParser.json())
app.use(bodyParser.urlencoded({ extended: true }))

// HTTP route for / to the index page.
app.get('/', function(req, res) {
    res.sendFile(path.join(__dirname, '/index.html'))
})

// HTTP route for getting a voronoi cell from vectors.
app.post('/voronoi_cell', function(req, res) {
    params = req.body
    // Running the algorithm script
    const options = {
        mode: 'text',
        pythonOptions: ['-u'], // get print results in real-time
        args: [params['x1'], params['y1'], params['z1'], params['x2'], params['y2'], params['z2'], params['x3'], params['y3'], params['z3'], params['m']]
    }
    PythonShell.run('algorithm.py', options, function (err, results) {
        if (err) 
            console.log(err)
        // Results is an array consisting of messages collected during execution
        res.send(results)
    })
})

// HTTP route for getting a voronoi cell from neigbors.
app.post('/voronoi_cell_by_points', function(req, res) {
    params = req.body
    const points = params['points'].split(' ')
    // Running the algorithm script
    const options = {
        mode: 'text',
        pythonOptions: ['-u'], // get print results in real-time
        args: ['points'].concat(points)
    }
    PythonShell.run('algorithm.py', options, function (err, results) {
        if (err) 
            console.log(err)
        // Results is an array consisting of messages collected during execution
        res.send(results)
    })   
})
