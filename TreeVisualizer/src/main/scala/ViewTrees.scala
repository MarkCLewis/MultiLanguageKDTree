import scalafx.Includes._
import scalafx.application.JFXApp3
import scalafx.scene.Scene
import scalafx.scene.paint.Color._
import scalafx.scene.shape.Rectangle
import scalafx.scene.canvas.Canvas
import scalafx.stage.FileChooser
import scalafx.scene.canvas.GraphicsContext

object ViewTrees extends JFXApp3 {
  val viewSize = 6.0
  override def start(): Unit = {
    stage = new JFXApp3.PrimaryStage {
      title.value = "Hello Stage"
      width = 1200
      height = 1200
      scene = new Scene {
        fill = White
        val canvas = new Canvas(1200, 1200)
        val gc = canvas.graphicsContext2D
        content = canvas
        canvas.onMouseClicked = e => {
          val chooser = new FileChooser
          val file = chooser.showOpenDialog(stage)
          println(file)
          if file != null then drawTree(gc, file)
        }
      }
    }
  }

  def drawTree(gc: GraphicsContext, file: java.io.File) = {
    val source = io.Source.fromFile(file)
    val lines = source.getLines()
    val num = lines.next().toInt
    gc.fill = White
    gc.fillRect(0, 0, gc.canvas.width(), gc.canvas.height())
    recur(gc, lines, -viewSize, viewSize, -viewSize, viewSize)
    source.close()
  }

  def recur(gc: GraphicsContext, lines: Iterator[String], minx: Double, maxx: Double, miny: Double, maxy: Double): Unit = {
    if lines.hasNext then {
      val line = lines.next()
      if line.startsWith("I") then {
        val Array(dim, value, left, right) = line.split(" +").drop(1).map(_.toDouble)
        if dim == 0 then {
          val px = (value + viewSize) / (2 * viewSize) * gc.canvas.width()
          val py1 = (miny + viewSize) / (2 * viewSize) * gc.canvas.height()
          val py2 = (maxy + viewSize) / (2 * viewSize) * gc.canvas.height()
          gc.stroke = Black
          gc.strokeLine(px, py1, px, py2)
          recur(gc, lines, minx, value, miny, maxy)
          recur(gc, lines, value, maxx, miny, maxy)
        } else if dim == 1 then {
          val px1 = (minx + viewSize) / (2 * viewSize) * gc.canvas.width()
          val px2 = (maxx + viewSize) / (2 * viewSize) * gc.canvas.width()
          val py = (value + viewSize) / (2 * viewSize) * gc.canvas.height()
          gc.stroke = Black
          gc.strokeLine(px1, py, px2, py)
          recur(gc, lines, minx, maxx, miny, value)
          recur(gc, lines, minx, maxx, value, maxy)
        } else {
          recur(gc, lines, minx, maxx, miny, maxy)
          recur(gc, lines, minx, maxx, miny, maxy)
        }
      } else if line.startsWith("L") then {
        val num = line.split(" +")(1).toInt
        for i <- 0 until num do {
          val Array(x, y, z) = lines.next().split(" +").map(_.toDouble)
          val px = (x + viewSize) / (2 * viewSize) * gc.canvas.width()
          val py = (y + viewSize) / (2 * viewSize) * gc.canvas.height()
          gc.stroke = Green
          gc.strokeLine(px, py, px, py)
        }
      } else {
        println("Bad file format!!!")
      }
    }
  }
}