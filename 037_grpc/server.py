from concurrent import futures

import converter_pb2
import converter_pb2_grpc
import grpc


class ConverterService(converter_pb2_grpc.ConverterServicer):
    def Convert(self, request, context):
        result = request.message.upper()
        return converter_pb2.ConvertResponse(result=result)


def serve():
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
    converter_pb2_grpc.add_ConverterServicer_to_server(ConverterService(), server)
    server.add_insecure_port("[::]:50051")
    server.start()
    print("Hello, gRPC!")
    server.wait_for_termination()


if __name__ == "__main__":
    serve()
